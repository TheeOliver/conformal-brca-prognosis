"""Uncertainty illustrations use hypothetical aggregates, never a selected patient."""

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose
from pandas.testing import assert_frame_equal

from brca.uncertainty import hypothetical_profiles, hypothetical_uncertainty_report
from brca.viz.uncertainty import export_hypothetical_uncertainty


def inputs():
    config = SimpleNamespace(
        feature_sets=SimpleNamespace(
            clinical=["age_at_diagnosis", "tumor_size", "tumor_grade"],
            clinical_molecular=["pam50_subtype"],
        ),
        preprocessing=SimpleNamespace(
            numeric_columns=["age_at_diagnosis", "tumor_size"], unknown_category="unknown"
        ),
        eda=SimpleNamespace(minimum_group_size=10, quantiles=[0.0, 0.25, 0.5, 0.75, 1.0]),
        bayesian=SimpleNamespace(
            ppc_time_grid_months=[12, 24, 60, 120, 300], predictive_interval_probability=0.9
        ),
        conformal=SimpleNamespace(alphas=[0.2, 0.1, 0.05]),
        outcome=SimpleNamespace(primary_endpoint="disease_specific"),
        figures=SimpleNamespace(
            width_inches=6.3, height_inches=4.0, font_size=9, line_width=1.6, dpi=80
        ),
    )
    train = pd.DataFrame(
        {
            "patient_id": [f"SYNTHETIC-{j}" for j in range(40)],
            "age_at_diagnosis": np.arange(40) + 35.0,
            "tumor_size": np.arange(40) + 10.0,
            "tumor_grade": ["3"] * 30 + ["2"] * 10,
            "pam50_subtype": ["LumA"] * 20 + ["Basal"] * 15 + ["small"] * 5,
            "event": [True] * 25 + [False] * 15,
            "survival_months": np.arange(40) + 1.0,
        }
    )
    return train, config


class PreparedProfiles:
    def transform(self, frame):
        assert "patient_id" not in frame and "event" not in frame
        return frame[["age_at_diagnosis", "tumor_size"]]


class AnalyticModel:
    """Known exponential distribution permits independent interval assertions."""

    def predict_survival_bands(self, X, times):
        scale = 2 * X["age_at_diagnosis"].to_numpy()[:, None]
        return {
            "mean": np.exp(-times / scale),
            "lower": np.exp(-times / (scale * 0.8)),
            "upper": np.exp(-times / (scale * 1.2)),
        }

    def predict_quantiles(self, X, probabilities):
        return -2 * X["age_at_diagnosis"].to_numpy()[:, None] * np.log1p(-probabilities)


def test_profiles_are_training_aggregates_exclude_small_levels_and_ignore_outcomes():
    train, config = inputs()
    profiles, definitions = hypothetical_profiles(train, config)
    assert len(profiles) == 4
    assert list(profiles.columns) == [*config.feature_sets.clinical, "pam50_subtype"]
    assert set(profiles.pam50_subtype) == {"LumA", "Basal"}
    assert all(row["hypothetical"] for row in definitions)
    assert profiles.loc["reference", "age_at_diagnosis"] == train.age_at_diagnosis.median()
    assert profiles.loc[
        "age_lower_quartile", "age_at_diagnosis"
    ] == train.age_at_diagnosis.quantile(0.25)
    assert profiles.loc[
        "age_upper_quartile", "age_at_diagnosis"
    ] == train.age_at_diagnosis.quantile(0.75)
    assert not set(profiles.index).intersection(train.patient_id)
    changed = train.assign(event=False, survival_months=9999)
    repeated, metadata = hypothetical_profiles(changed, config)
    assert_frame_equal(profiles, repeated)
    assert definitions == metadata


def test_predictive_time_intervals_follow_analytic_distribution_and_remain_unrestricted(tmp_path):
    train, config = inputs()
    profiles, definitions = hypothetical_profiles(train, config)
    models = {key: AnalyticModel() for key in ("clinical", "clinical_molecular")}
    preprocessors = {key: PreparedProfiles() for key in models}
    report = hypothetical_uncertainty_report(profiles, definitions, models, preprocessors, config)
    assert report["patient_records_included"] is False
    assert report["test_accessed"] is False
    assert report["calibration_accessed"] is False
    assert report["credible_band_type"] == "pointwise_parameter_uncertainty_in_survival_function"
    reference = report["profiles"][0]["models"]["clinical"]
    scale = 2 * profiles.loc["reference", "age_at_diagnosis"]
    assert_allclose(reference["predictive_median_months"], scale * np.log(2))
    for interval in reference["predictive_intervals"]:
        assert_allclose(interval["lower_months"], -scale * np.log1p(-interval["alpha"] / 2))
        assert_allclose(interval["upper_months"], -scale * np.log(interval["alpha"] / 2))
        assert (
            interval["lower_months"]
            < reference["predictive_median_months"]
            < interval["upper_months"]
        )
    assert reference["predictive_intervals"][-1]["upper_months"] > max(report["time_grid_months"])
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps(report, allow_nan=False))
    artifacts = export_hypothetical_uncertainty(path, tmp_path / "figures", config)
    assert len(artifacts) == 2 * len(profiles)
    assert all(artifact.stat().st_size > 1000 for artifact in artifacts)
