"""Descriptive exploration restricted to the frozen training cohort."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.util import Surv


def outcome_labels(cohort: pd.DataFrame, endpoint: str) -> np.ndarray:
    column = f"event_{endpoint}"
    if column not in cohort:
        raise ValueError("The cohort lacks the requested registered endpoint.")
    return Surv.from_arrays(
        cohort[column].to_numpy(bool), cohort["survival_months"].to_numpy(float)
    )


def _km_summary(frame: pd.DataFrame, endpoint: str, config) -> dict:
    y = outcome_labels(frame, endpoint)
    time, survival, interval = kaplan_meier_estimator(
        y["event"], y["time"], conf_type="log-log", conf_level=config.eda.km_confidence_level
    )
    followup_time, followup_survival = kaplan_meier_estimator(y["event"], y["time"], reverse=True)
    median_indices = np.flatnonzero(survival <= 0.5)
    followup_indices = np.flatnonzero(followup_survival <= 0.5)
    return {
        "n": len(frame),
        "events": int(y["event"].sum()),
        "censored": int((~y["event"]).sum()),
        "event_fraction": float(y["event"].mean()),
        "median_survival_months": float(time[median_indices[0]]) if len(median_indices) else None,
        "reverse_km_median_followup_months": float(followup_time[followup_indices[0]])
        if len(followup_indices)
        else None,
        "times_months": time.tolist(),
        "survival": survival.tolist(),
        "lower": interval[0].tolist(),
        "upper": interval[1].tolist(),
        "at_risk": {
            str(t): int((y["time"] >= t).sum()) for t in config.eda.followup_horizons_months
        },
        "censor_times_months": np.unique(y["time"][~y["event"]]).tolist(),
    }


def training_eda(cohort: pd.DataFrame, config) -> dict:
    """Return aggregate descriptions, never transformed or individual patient rows."""
    report = {
        "split": "train",
        "n_training": len(cohort),
        "numeric": {},
        "categorical": {},
        "survival": {},
    }
    for column in config.preprocessing.numeric_columns:
        values = cohort[column].dropna().to_numpy(float)
        counts, edges = np.histogram(values, bins=config.eda.histogram_bins)
        report["numeric"][column] = {
            "n_observed": len(values),
            "n_missing": int(cohort[column].isna().sum()),
            "quantiles": {str(q): float(np.quantile(values, q)) for q in config.eda.quantiles},
            "histogram_edges": edges.tolist(),
            "histogram_counts": counts.tolist(),
        }
    for column in config.preprocessing.categorical_columns:
        counts = (
            cohort[column]
            .astype(object)
            .fillna(config.preprocessing.unknown_category)
            .astype(str)
            .value_counts()
        )
        report["categorical"][column] = {
            str(level): {"n": int(n), "small_group": bool(n < config.eda.minimum_group_size)}
            for level, n in counts.items()
        }
    for endpoint in [config.outcome.primary_endpoint, *config.outcome.sensitivity_endpoints]:
        summary = {"overall": _km_summary(cohort, endpoint, config), "groups": {}}
        for column in config.conformal.coverage_subgroups:
            summary["groups"][column] = {}
            for level, group in cohort.groupby(column, dropna=False, observed=True):
                summary["groups"][column][str(level)] = (
                    _km_summary(group, endpoint, config)
                    if len(group) >= config.eda.minimum_group_size
                    else {"n": len(group), "status": "suppressed_small_group"}
                )
        report["survival"][endpoint] = summary
    report["numeric_spearman_correlation"] = (
        cohort[config.preprocessing.numeric_columns].corr(method="spearman").to_dict()
    )
    report["missingness_mechanism"] = (
        "These descriptions cannot establish that missingness is random."
    )
    return report
