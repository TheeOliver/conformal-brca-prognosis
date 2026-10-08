"""Explicit clinical source mapping and cohort exclusions; no fitted transforms."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

from brca.data.schema import validate

NUMERIC_MAPPING = {
    "age_at_diagnosis": "AGE_AT_DIAGNOSIS",
    "tumor_size": "TUMOR_SIZE",
    "tumor_grade": "GRADE",
    "lymph_nodes_positive": "LYMPH_NODES_EXAMINED_POSITIVE",
}
RECEPTOR_MAPPING = {
    "er_status": "ER_STATUS",
    "pr_status": "PR_STATUS",
    "her2_status": "HER2_STATUS",
}
SUBTYPES = {"LumA", "LumB", "Her2", "Basal", "Normal", "claudin-low", "NC"}
VITAL_STATES = {"Living", "Died of Disease", "Died of Other Causes"}
OS_STATES = {"0:LIVING", "1:DECEASED"}


def numeric_source(column: pd.Series, name: str) -> pd.Series:
    values = pd.to_numeric(column, errors="coerce")
    if (column.notna() & values.isna()).any() or np.isinf(values).any():
        raise ValueError(f"Invalid non-numeric or infinite source values in {name}.")
    return values.astype(float)


def prepare_cohort(raw: pd.DataFrame, config: SimpleNamespace) -> tuple[pd.DataFrame, dict]:
    """Define the DSS cohort once; retain the same patients for OS sensitivity."""
    required = (
        set(NUMERIC_MAPPING.values())
        | set(RECEPTOR_MAPPING.values())
        | {"OS_MONTHS", "OS_STATUS", "VITAL_STATUS", "CLAUDIN_SUBTYPE"}
    )
    if not required <= set(raw.columns):
        raise ValueError("The clinical snapshot lacks required source attributes.")
    if raw.index.has_duplicates or raw.index.isna().any():
        raise ValueError("Patient identifiers must be present and unique.")
    if config.outcome.primary_endpoint != "disease_specific":
        raise ValueError("The registered primary cohort uses disease-specific survival.")
    for field, levels in (
        ("OS_STATUS", OS_STATES),
        ("VITAL_STATUS", VITAL_STATES),
        ("CLAUDIN_SUBTYPE", SUBTYPES),
    ):
        if not raw[field].dropna().isin(levels).all():
            raise ValueError(f"Unexpected source category in {field}; mapping requires review.")
    time = numeric_source(raw["OS_MONTHS"], "OS_MONTHS")
    keep = pd.Series(True, index=raw.index)
    exclusions = {}
    for reason, invalid in (
        ("missing_overall_outcome", time.isna() | raw["OS_STATUS"].isna()),
        ("missing_cause_of_death", raw["VITAL_STATUS"].isna()),
        ("nonpositive_followup", time.le(0)),
    ):
        exclusions[reason] = int((keep & invalid).sum())
        keep &= ~invalid
    source = raw.loc[keep]
    inconsistent = source["OS_STATUS"].eq("0:LIVING") != source["VITAL_STATUS"].eq("Living")
    if inconsistent.any():
        raise ValueError("Overall status and recorded cause of death are inconsistent.")
    cohort = pd.DataFrame(index=source.index)
    cohort["patient_id"] = source.index.astype(str)
    cohort["survival_months"] = time.loc[keep]
    cohort["event_disease_specific"] = source["VITAL_STATUS"].eq("Died of Disease")
    cohort["event_overall"] = source["OS_STATUS"].eq("1:DECEASED")
    horizon = config.outcome.administrative_censoring_months
    if horizon is not None:
        if not np.isfinite(horizon) or horizon <= 0:
            raise ValueError("Administrative horizon must be positive or null.")
        before_horizon = cohort["survival_months"].le(horizon)
        for column in ("event_disease_specific", "event_overall"):
            cohort[column] &= before_horizon
        cohort["survival_months"] = cohort["survival_months"].clip(upper=horizon)
    cohort["event"] = cohort["event_disease_specific"]
    for name, field in NUMERIC_MAPPING.items():
        cohort[name] = numeric_source(source[field], field)
    invalid_predictors = (
        cohort["age_at_diagnosis"].dropna().le(0).any()
        or cohort["tumor_size"].dropna().lt(0).any()
        or cohort["lymph_nodes_positive"].dropna().lt(0).any()
        or not cohort["tumor_grade"].dropna().isin([1, 2, 3]).all()
    )
    if invalid_predictors:
        raise ValueError("Impossible clinical predictor coding; review the source mapping.")
    for name, field in RECEPTOR_MAPPING.items():
        if not source[field].dropna().isin(["Positive", "Negative"]).all():
            raise ValueError(f"Unexpected source category in {field}.")
        cohort[name] = source[field].str.lower()
    if not config.preprocessing.retain_claudin_low:
        raise ValueError("The registered feature comparison retains claudin-low.")
    cohort["pam50_subtype"] = (
        source["CLAUDIN_SUBTYPE"]
        .replace("NC", config.preprocessing.unknown_category)
        .fillna(config.preprocessing.unknown_category)
    )
    cohort = cohort.sort_index().reset_index(drop=True)
    validate(cohort)
    report = {
        "input_kind": "real_metabric",
        "raw_patients": len(raw),
        "exclusions_sequential": exclusions,
        "analysis_patients": len(cohort),
        "disease_specific_events": int(cohort["event_disease_specific"].sum()),
        "overall_events": int(cohort["event_overall"].sum()),
        "administrative_censoring_months": horizon,
        "primary_endpoint": config.outcome.primary_endpoint,
        "sensitivity_endpoints": config.outcome.sensitivity_endpoints,
        "missing_predictors": {
            name: int(cohort[name].isna().sum()) for name in (*NUMERIC_MAPPING, *RECEPTOR_MAPPING)
        },
        "subtype_counts": cohort["pam50_subtype"].value_counts().to_dict(),
        "fit_transforms": False,
    }
    return cohort, report


def raw_audit(raw: pd.DataFrame, config: SimpleNamespace) -> dict:
    """Aggregate source audit before splitting; no associations or model selection."""
    from copy import deepcopy

    audit_config = deepcopy(config)
    audit_config.outcome.administrative_censoring_months = None
    cohort, report = prepare_cohort(raw, audit_config)
    times = cohort["survival_months"]
    report["observed_followup_quantiles_months"] = {
        str(q): float(times.quantile(q)) for q in config.eda.quantiles
    }
    report["followup_tail"] = {
        str(horizon): {
            "observed_beyond": int(times.gt(horizon).sum()),
            "dss_events_beyond": int((times.gt(horizon) & cohort["event_disease_specific"]).sum()),
            "os_events_beyond": int((times.gt(horizon) & cohort["event_overall"]).sum()),
        }
        for horizon in config.eda.followup_horizons_months
    }
    predictors = list(config.feature_sets.clinical)
    report["any_missing_predictor"] = int(cohort[predictors].isna().any(axis=1).sum())
    report["complete_cases"] = int(cohort[predictors].notna().all(axis=1).sum())
    report["scope"] = "aggregate source/cohort audit before freezing the real split"
    return report
