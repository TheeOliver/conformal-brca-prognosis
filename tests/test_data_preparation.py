"""Cohort mapping and preprocessing on fabricated source tables only."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from brca.config import load_config
from brca.data.cbioportal import save_raw
from brca.data.load import load_cohort, load_raw_clinical
from brca.data.prepare import prepare_cohort, raw_audit
from brca.data.preprocessing import TrainingPreprocessor

FIXTURE = Path(__file__).parent / "fixtures" / "synthetic_metabric.csv"


@pytest.fixture
def raw():
    n = 20
    return pd.DataFrame(
        {
            "OS_MONTHS": [str(12 + i * 20) for i in range(n)],
            "OS_STATUS": ["1:DECEASED", "0:LIVING"] * (n // 2),
            "VITAL_STATUS": ["Died of Disease", "Living", "Died of Other Causes", "Living"]
            * (n // 4),
            "AGE_AT_DIAGNOSIS": ["55"] * n,
            "TUMOR_SIZE": ["20"] * n,
            "GRADE": ["2"] * n,
            "LYMPH_NODES_EXAMINED_POSITIVE": ["0"] * n,
            "ER_STATUS": ["Positive"] * n,
            "PR_STATUS": ["Negative"] * n,
            "HER2_STATUS": ["Negative"] * n,
            "CLAUDIN_SUBTYPE": ["LumA", "LumB", "claudin-low", "NC"] * (n // 4),
        },
        index=pd.Index([f"FAKE-{i:03d}" for i in range(n)], name="patient_id"),
    )


def test_mapping_retains_competing_deaths_and_both_endpoint_labels(raw):
    cohort, report = prepare_cohort(raw, load_config())
    assert len(cohort) == len(raw)
    assert cohort.event.dtype == bool
    assert cohort.event.sum() == 5
    assert cohort.event_overall.sum() == 10
    assert cohort.event.equals(cohort.event_disease_specific)
    assert set(cohort.pam50_subtype) == {"LumA", "LumB", "claudin-low", "unknown"}
    assert report["fit_transforms"] is False
    assert cohort.survival_months.max() > 300


def test_exclusions_are_sequential_and_missing_predictors_retained(raw):
    raw.loc[raw.index[0], "OS_MONTHS"] = None
    raw.loc[raw.index[1], "VITAL_STATUS"] = None
    raw.loc[raw.index[2], "OS_MONTHS"] = "0"
    raw.loc[raw.index[3], "TUMOR_SIZE"] = None
    cohort, report = prepare_cohort(raw, load_config())
    assert len(cohort) == len(raw) - 3
    assert report["exclusions_sequential"] == {
        "missing_overall_outcome": 1,
        "missing_cause_of_death": 1,
        "nonpositive_followup": 1,
    }
    assert cohort.tumor_size.isna().sum() == 1


def test_admin_censoring_is_explicit_and_correct_for_both_endpoints(raw):
    config = load_config()
    config.outcome.administrative_censoring_months = 100
    cohort, report = prepare_cohort(raw, config)
    assert cohort.survival_months.max() == 100
    assert cohort.event.sum() == 2
    assert cohort.event_overall.sum() == 3
    assert report["administrative_censoring_months"] == 100


@pytest.mark.parametrize(
    "field,value",
    [
        ("OS_STATUS", "unexpected"),
        ("GRADE", "9"),
        ("OS_MONTHS", "bad"),
        ("TUMOR_SIZE", "-1"),
        ("CLAUDIN_SUBTYPE", "unmapped"),
    ],
)
def test_invalid_source_coding_fails_without_printing_patients(raw, field, value):
    raw.loc[raw.index[0], field] = value
    with pytest.raises(ValueError) as error:
        prepare_cohort(raw, load_config())
    assert "FAKE-" not in str(error.value)


def test_raw_audit_is_aggregate_and_does_not_apply_admin_censoring(raw):
    config = load_config()
    config.outcome.administrative_censoring_months = 100
    audit = raw_audit(raw, config)
    assert audit["observed_followup_quantiles_months"]["1.0"] == 392
    assert audit["followup_tail"]["300"]["observed_beyond"] == 5
    assert "FAKE-" not in str(audit)


def test_checksummed_raw_loader_joins_on_patient_id(tmp_path):
    payloads = {
        "study": {"studyId": "synthetic", "publicStudy": True},
        "attributes": [],
        "patient": [
            {"patientId": "FAKE-1", "clinicalAttributeId": "AGE_AT_DIAGNOSIS", "value": "55"}
        ],
        "sample": [
            {
                "patientId": "FAKE-1",
                "sampleId": "DIFFERENT-SAMPLE-ID",
                "clinicalAttributeId": "GRADE",
                "value": "2",
            }
        ],
    }
    save_raw(payloads, tmp_path)
    raw = load_raw_clinical(tmp_path)
    assert len(raw) == 1 and set(raw.columns) == {"AGE_AT_DIAGNOSIS", "GRADE"}
    source = tmp_path / "clinical_patient.json.gz"
    source.write_bytes(source.read_bytes() + b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        load_raw_clinical(tmp_path)


def test_preprocessing_learns_only_training_values_and_fixed_categories():
    cohort = load_cohort(FIXTURE)
    train = cohort.iloc[:100].copy()
    test = cohort.iloc[100:110].copy()
    config = load_config()
    preprocessor = TrainingPreprocessor(config, "clinical_molecular")
    preprocessor.fit(train)
    before = deepcopy((preprocessor.medians_, preprocessor.means_, preprocessor.categories_))
    test.loc[test.index[0], "tumor_size"] = np.nan
    test.loc[test.index[1], "pam50_subtype"] = "unseen-synthetic-category"
    test.loc[test.index[2], "age_at_diagnosis"] = 500
    output = preprocessor.transform(test)
    assert before == (preprocessor.medians_, preprocessor.means_, preprocessor.categories_)
    assert np.isfinite(output.to_numpy()).all()
    assert output.loc[test.index[0], "tumor_size_missing"] == 1
    assert output.loc[test.index[1], "pam50_subtype=unknown"] == 1
    assert list(output) == list(preprocessor.transform(train))


def test_clinical_preprocessing_is_identical_between_feature_arms():
    cohort = load_cohort(FIXTURE)
    train = cohort.iloc[:100]
    clinical = TrainingPreprocessor(load_config(), "clinical").fit_transform(train)
    molecular = TrainingPreprocessor(load_config(), "clinical_molecular").fit_transform(train)
    pd.testing.assert_frame_equal(clinical, molecular[clinical.columns])


def test_preprocessor_roundtrip(tmp_path):
    cohort = load_cohort(FIXTURE)
    preprocessor = TrainingPreprocessor(load_config(), "clinical").fit(cohort)
    path = tmp_path / "transform.pkl"
    preprocessor.save(path)
    pd.testing.assert_frame_equal(
        preprocessor.transform(cohort), TrainingPreprocessor.load(path).transform(cohort)
    )
