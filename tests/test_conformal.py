"""Synthetic coverage contracts; these numbers are validation, never thesis results."""

import json

import numpy as np
import pytest
from sksurv.util import Surv

from brca.conformal.survival import (
    CensoringDistribution,
    ConformalCalibration,
    PredictionIntervals,
    _draw_summary,
    calibrate_intervals,
    calibrate_ipcw_intervals,
    finite_sample_quantile,
    paired_interval_comparison,
    score_bounds,
    subgroup_coverage,
)
from brca.conformal.survival import (
    coverage_summary as _coverage_summary,
)


def coverage_summary(*args, **kwargs):
    kwargs.setdefault("minimum_bootstrap_success_fraction", 0.9)
    return _coverage_summary(*args, **kwargs, confidence_level=0.95, min_censoring_survival=0.05)


def labels(time, event=None):
    time = np.asarray(time, float)
    event = np.ones(len(time), bool) if event is None else np.asarray(event, bool)
    return Surv.from_arrays(event=event, time=time)


@pytest.mark.parametrize("alpha,k", [(0.2, 80), (0.1, 90), (0.05, 95)])
def test_exact_order_statistic(alpha, k):
    scores = np.arange(1, 100, dtype=float)
    assert finite_sample_quantile(scores, alpha) == k
    assert finite_sample_quantile(scores[::-1], alpha) == k


def test_small_and_empty_calibration_require_infinity():
    assert np.isinf(finite_sample_quantile(np.arange(1, 10), 0.05))
    assert np.isinf(finite_sample_quantile(np.array([]), 0.05))
    assert finite_sample_quantile(np.arange(1, 20), 0.05) == 19


@pytest.mark.parametrize("alpha", [0, 1, -0.1, np.nan])
def test_invalid_alpha(alpha):
    with pytest.raises(ValueError, match="alpha"):
        finite_sample_quantile(np.array([1]), alpha)


def test_observed_and_censored_score_envelopes():
    y = labels([20, 20, 110, 50], [True, False, False, True])
    scores = score_bounds(y, np.full(4, 30), np.full(4, 70), horizon=100)
    np.testing.assert_array_equal(scores, [10, 30, 30, 0])
    assert np.isinf(score_bounds(y, np.full(4, 30), np.full(4, 70), horizon=np.inf)[1])


@pytest.mark.parametrize("alpha", [0.2, 0.1, 0.05])
def test_exchangeable_uncensored_coverage_across_repeated_calibrations(alpha):
    rng = np.random.default_rng(39142)
    n_cal, n_test, repeats = 199, 300, 250
    rates = []
    for _ in range(repeats):
        latent = rng.uniform(1, 100, n_cal + n_test)
        calibration = calibrate_intervals(
            labels(latent[:n_cal]), np.full(n_cal, 50), np.full(n_cal, 50), alpha, horizon=100
        )
        intervals = calibration.predict(np.full(n_test, 50), np.full(n_test, 50))
        rates.append(
            np.mean((intervals.lower <= latent[n_cal:]) & (latent[n_cal:] <= intervals.upper))
        )
    # The Monte Carlo unit is the independent calibration/test replicate, not
    # each correlated test patient sharing a calibrated quantile.
    expected = np.ceil((n_cal + 1) * (1 - alpha)) / (n_cal + 1)
    mc_se = np.std(rates, ddof=1) / np.sqrt(repeats)
    assert abs(np.mean(rates) - expected) < 4 * mc_se + 0.001


@pytest.mark.parametrize("alpha", [0.2, 0.1, 0.05])
def test_dependent_censoring_and_stratification_preserve_conservative_coverage(alpha):
    rng = np.random.default_rng(394)
    rates = []
    for _ in range(100):
        latent = rng.uniform(1, 150, 1400)
        # Deliberately dependent censoring: worst for long survival.
        censor = np.where(latent > 60, rng.uniform(15, 110, len(latent)), 170)
        event = latent <= censor
        y = labels(np.minimum(latent, censor), event)
        cal_idx, test_idx = [], []
        for level in [False, True]:
            idx = np.flatnonzero(event == level)
            rng.shuffle(idx)
            cal_idx.extend(idx[:150])
            test_idx.extend(idx[150:350])
        cal_idx, test_idx = np.array(cal_idx), np.array(test_idx)
        calibration = calibrate_intervals(
            y[cal_idx],
            np.full(len(cal_idx), 25),
            np.full(len(cal_idx), 60),
            alpha,
            horizon=100,
            strata=event[cal_idx],
            expected_strata=[False, True],
        )
        intervals = calibration.predict(np.full(len(test_idx), 25), np.full(len(test_idx), 60))
        target = np.minimum(latent[test_idx], 100)
        rates.append(np.mean((intervals.lower <= target) & (target <= intervals.upper)))
    mc_se = np.std(rates, ddof=1) / np.sqrt(len(rates))
    assert np.mean(rates) >= 1 - alpha - 4 * mc_se


def test_heavy_censoring_unrestricted_infinity_and_restricted_full_support():
    y = labels(np.full(100, 1), np.zeros(100, bool))
    lower, upper = np.full(100, 30), np.full(100, 40)
    state = calibrate_intervals(y, lower, upper, 0.1, horizon=np.inf)
    intervals = state.predict(lower[:3], upper[:3])
    assert np.all(np.isinf(intervals.width))
    assert np.all(intervals.lower == 0)
    restricted = calibrate_intervals(y, lower, upper, 0.1, horizon=100).predict(lower, upper)
    assert np.all(restricted.lower == 0)
    assert np.all(restricted.upper == 100)


def test_missing_stratum_returns_full_support_and_serializes():
    y = labels(np.arange(1, 51))
    state = calibrate_intervals(
        y,
        np.zeros(50),
        np.ones(50),
        0.1,
        horizon=100,
        strata=np.ones(50, bool),
        expected_strata=[False, True],
    )
    assert state.stratum_counts["False"] == 0
    assert np.isinf(state.quantile)
    restored = ConformalCalibration.from_dict(
        json.loads(json.dumps(state.to_dict(), allow_nan=False))
    )
    assert np.isinf(restored.quantile)
    assert restored.predict(np.array([20]), np.array([30])).width[0] == 100


def test_reverse_km_event_censor_ties_and_support():
    censoring = CensoringDistribution.fit(labels([1, 1, 2, 3], [True, False, False, True]))
    # At time 1 the event is removed before a censoring jump: 1 - 1/3.
    np.testing.assert_allclose(
        censoring.survival_ge(np.array([1, 1.5, 2, 2.5])), [1, 2 / 3, 2 / 3, 1 / 3]
    )
    assert np.isnan(censoring.survival_ge(np.array([4]))[0])


def test_coverage_bounds_cover_realised_latent_coverage():
    latent = np.array([20, 80, 140, 50, 100], float)
    y = labels([20, 40, 90, 50, 70], [True, False, False, True, False])
    intervals = PredictionIntervals(
        np.array([10, 30, 95, 55, 0]),
        np.array([25, 100, 100, 80, 60]),
        alpha=0.1,
        n_cal=100,
        horizon=100,
    )
    report = coverage_summary(y, intervals, labels([150, 160]), bootstrap_replicates=50, seed=3)
    target = np.minimum(latent, 100)
    actual = np.mean((intervals.lower <= target) & (target <= intervals.upper))
    assert report["observable_coverage_lower"] == 0.4
    assert report["observable_coverage_upper"] == 0.6
    assert report["observable_coverage_lower"] <= actual <= report["observable_coverage_upper"]
    assert report["width_median_months"] == 25
    json.dumps(report, allow_nan=False)


def test_censoring_at_upper_endpoint_is_not_possibly_covered():
    y = labels([50], [False])
    intervals = PredictionIntervals(np.array([0]), np.array([50]), 0.1, 100, 100)
    report = coverage_summary(y, intervals, labels([150]), bootstrap_replicates=0, seed=3)
    assert report["observable_coverage_upper"] == 0
    assert report["ipcw_status"] == "unavailable_no_observed_targets"


def test_known_restricted_horizon_counts_censored_followup():
    y = labels([100, 130, 50], [False, False, True])
    intervals = PredictionIntervals(np.array([90, 90, 45]), np.array([100, 100, 55]), 0.1, 100, 100)
    report = coverage_summary(y, intervals, labels([150, 160]), bootstrap_replicates=30, seed=3)
    assert report["n_known_target"] == 3
    assert report["observable_coverage_lower"] == report["observable_coverage_upper"] == 1
    assert report["ipcw_coverage"] == 1


def test_ipcw_estimate_recovers_known_coverage_under_independent_censoring():
    rng = np.random.default_rng(9342)
    n = 18000
    latent_train, latent_test = rng.uniform(1, 120, (2, n))
    censor_train, censor_test = rng.uniform(1, 180, (2, n))
    y_train = labels(np.minimum(latent_train, censor_train), latent_train <= censor_train)
    y_test = labels(np.minimum(latent_test, censor_test), latent_test <= censor_test)
    intervals = PredictionIntervals(np.full(n, 12), np.full(n, 95), 0.2, 400, 100)
    report = coverage_summary(y_test, intervals, y_train, bootstrap_replicates=50, seed=23)
    target = np.minimum(latent_test, 100)
    truth = np.mean((intervals.lower <= target) & (target <= intervals.upper))
    assert abs(report["ipcw_coverage"] - truth) < 0.025
    assert abs(report["ipcw_ht_coverage"] - truth) < 0.025
    assert report["ipcw_coverage_ci"][0] <= report["ipcw_coverage"] <= report["ipcw_coverage_ci"][1]


@pytest.mark.parametrize("alpha", [0.2, 0.1, 0.05])
def test_ipcw_calibration_reduces_to_exact_quantile_without_censoring(alpha):
    y = labels(np.arange(1, 100))
    lo, hi = np.full(99, 50), np.full(99, 50)
    expected = calibrate_intervals(y, lo, hi, alpha, horizon=100)
    actual = calibrate_ipcw_intervals(
        y, lo, hi, labels([150, 160]), alpha, horizon=100, min_censoring_survival=0.05
    )
    assert expected.quantile == actual.quantile


def test_ipcw_calibration_censored_coverage_simulation():
    rng = np.random.default_rng(9805)
    rates = []
    for _ in range(150):
        t = rng.uniform(1, 120, 2200)
        c = rng.uniform(1, 180, 2200)
        y = labels(np.minimum(t, c), t <= c)
        state = calibrate_ipcw_intervals(
            y[1000:1600],
            np.full(600, 50),
            np.full(600, 50),
            y[:1000],
            0.1,
            horizon=115,
            min_censoring_survival=0.05,
        )
        intervals = state.predict(np.full(600, 50), np.full(600, 50))
        z = np.minimum(t[1600:], 115)
        rates.append(np.mean((intervals.lower <= z) & (z <= intervals.upper)))
    # Asymptotic estimated-weight sensitivity, tested under its own assumptions.
    assert 0.87 <= np.mean(rates) <= 0.96


def test_missing_censoring_support_does_not_clip_or_extrapolate():
    y_train = labels([10, 20, 30], [False, False, False])
    y_test = labels([15, 25], [True, False])
    intervals = PredictionIntervals(np.zeros(2), np.full(2, 100), 0.1, 100, 100)
    report = coverage_summary(y_test, intervals, y_train, bootstrap_replicates=10, seed=3)
    assert report["ipcw_status"] == "unavailable_censoring_support"
    assert report["ipcw_coverage"] is None
    assert report["observable_coverage_lower"] == 1
    assert report["observable_coverage_lower_ci"] == [1, 1]
    assert report["width_median_ci_months"] == [100, 100]
    assert report["full_support_fraction_ci"] == [1, 1]
    assert report["infinite_interval_fraction_ci"] == [0, 0]
    state = calibrate_ipcw_intervals(
        y_test, np.zeros(2), np.ones(2), y_train, 0.1, horizon=100, min_censoring_survival=0.05
    )
    assert np.isinf(state.quantile)


def test_missing_subgroup_is_reported_with_null_estimates():
    y = labels([10, 20])
    intervals = PredictionIntervals(np.zeros(2), np.full(2, 100), 0.1, 100, 100)
    report = subgroup_coverage(
        y,
        intervals,
        labels([150, 160]),
        groups={"pam50": np.array(["A", "A"])},
        expected_levels={"pam50": ["A", "B"]},
        bootstrap_replicates=10,
        seed=1,
        confidence_level=0.95,
        min_censoring_survival=0.05,
        minimum_bootstrap_success_fraction=0.9,
    )
    missing = report["pam50"]["B"]
    assert missing["n_test"] == 0
    assert missing["observable_coverage_lower"] is None
    assert missing["ipcw_status"] == "unavailable_empty_subgroup"
    json.dumps(report, allow_nan=False)


def test_infinite_widths_serialize_explicitly_without_nan():
    intervals = PredictionIntervals(np.zeros(2), np.full(2, np.inf), 0.05, 2, np.inf)
    report = coverage_summary(
        labels([1, 2]), intervals, labels([10]), bootstrap_replicates=0, seed=1
    )
    assert report["infinite_interval_fraction"] == 1
    assert report["width_median_is_infinite"]
    assert report["width_median_months"] is None
    json.dumps(report, allow_nan=False)


def test_width_bootstrap_retains_infinity_flags_when_ipcw_is_unavailable():
    intervals = PredictionIntervals(np.zeros(2), np.full(2, np.inf), 0.05, 2, np.inf)
    report = coverage_summary(
        labels([1, 2]), intervals, labels([10]), bootstrap_replicates=20, seed=1
    )
    assert report["width_median_ci_months"] == [None, None]
    assert report["bootstrap_statistics"]["width_median_months"]["ci_is_infinite"] == [True, True]
    assert report["infinite_interval_fraction_ci"] == [1, 1]
    assert report["ipcw_coverage_ci"] is None
    json.dumps(report, allow_nan=False)


def test_paired_width_difference_uses_the_identical_patient_draws():
    y = labels([10, 20, 30, 40])
    clinical = PredictionIntervals(np.zeros(4), np.array([15, 20, 25, 40]), 0.1, 100, 100)
    molecular = PredictionIntervals(np.zeros(4), clinical.upper + 10, 0.1, 100, 100)
    report = paired_interval_comparison(
        y,
        clinical,
        molecular,
        labels([150, 160]),
        bootstrap_replicates=100,
        seed=3,
        confidence_level=0.95,
        min_censoring_survival=0.05,
        minimum_bootstrap_success_fraction=0.9,
    )
    differences = report["differences"]
    assert differences["width_median_months"]["estimate"] == 10
    # A constant within-person difference has zero paired bootstrap variation.
    assert differences["width_median_months"]["ci"] == [10, 10]
    assert differences["ipcw_coverage"]["estimate"] == 0.25
    assert differences["ipcw_coverage"]["valid_replicates"] == 100
    assert differences["observable_coverage_lower"]["estimate"] == 0.25
    json.dumps(report, allow_nan=False)


def test_paired_comparison_still_reports_width_without_censoring_support():
    y = labels([10, 20])
    clinical = PredictionIntervals(np.zeros(2), np.full(2, 30), 0.1, 100, 100)
    molecular = PredictionIntervals(np.zeros(2), np.full(2, 40), 0.1, 100, 100)
    report = paired_interval_comparison(
        y,
        clinical,
        molecular,
        labels([10, 15]),
        bootstrap_replicates=20,
        seed=3,
        confidence_level=0.95,
        min_censoring_survival=0.05,
        minimum_bootstrap_success_fraction=0.9,
    )
    assert report["ipcw_status"] == "unavailable_censoring_support"
    assert report["differences"]["ipcw_coverage"]["estimate"] is None
    assert report["differences"]["width_median_months"]["ci"] == [10, 10]


def test_paired_comparison_rejects_different_estimands():
    y = labels([10, 20])
    clinical = PredictionIntervals(np.zeros(2), np.full(2, 30), 0.1, 100, 120)
    molecular = PredictionIntervals(np.zeros(2), np.full(2, 40), 0.1, 100, 300)
    with pytest.raises(ValueError, match="same alpha and horizon"):
        paired_interval_comparison(
            y,
            clinical,
            molecular,
            labels([150]),
            bootstrap_replicates=20,
            seed=3,
            confidence_level=0.95,
            min_censoring_survival=0.05,
            minimum_bootstrap_success_fraction=0.9,
        )


def test_horizons_are_passed_explicitly_and_not_selected_from_support():
    y = labels([150, 150], [False, False])
    lower, upper = np.full(2, 25), np.full(2, 100)
    np.testing.assert_array_equal(score_bounds(y, lower, upper, horizon=120), [20, 20])
    np.testing.assert_array_equal(score_bounds(y, lower, upper, horizon=300), [200, 200])


def test_bootstrap_requires_two_valid_draws_and_flags_partial_success():
    one = _draw_summary(np.array([0.5, np.nan]), 0.95, 0.4)
    assert one["valid_replicates"] == 1
    assert one["failed_replicates"] == 1
    assert one["status"] == "bootstrap_gate_failed"
    assert one["ci"] is None
    partial = _draw_summary(np.r_[np.arange(9), np.nan], 0.95, 0.9)
    assert partial["status"] == "conditional_on_successful_resamples"
    assert partial["approximate_due_to_failed_resamples"]
    assert partial["valid_replicates"] == 9 and partial["failed_replicates"] == 1
    assert partial["ci"] is not None
    failed = _draw_summary(np.r_[np.arange(8), np.nan, np.nan], 0.95, 0.9)
    assert failed["status"] == "bootstrap_gate_failed"
    assert failed["ci"] is None


def test_bootstrap_preserves_primary_event_counts_for_secondary_endpoint():
    # All events here are OS deaths, while only the first two were DSS events.
    y = labels([10, 20, 30, 40])
    primary_strata = np.array([True, True, False, False])
    intervals = PredictionIntervals(np.zeros(4), np.array([15, 25, 25, 30]), 0.1, 100, 100)
    report = coverage_summary(
        y,
        intervals,
        labels([150]),
        bootstrap_replicates=100,
        seed=3,
        bootstrap_strata=primary_strata,
    )
    assert report["bootstrap_sampling"] == "within_primary_event_strata"
    assert report["observable_coverage_lower_ci"] == [0.5, 0.5]
    assert report["ipcw_coverage_ci"] == [0.5, 0.5]


def test_subgroup_bootstrap_subsets_primary_event_strata():
    y = labels([10, 20, 30, 40])
    intervals = PredictionIntervals(np.zeros(4), np.array([15, 15, 35, 35]), 0.1, 100, 100)
    report = subgroup_coverage(
        y,
        intervals,
        labels([150]),
        groups={"subtype": np.array(["A", "A", "B", "B"])},
        bootstrap_strata=np.array([True, False, True, False]),
        bootstrap_replicates=30,
        seed=9,
        confidence_level=0.95,
        min_censoring_survival=0.05,
        minimum_bootstrap_success_fraction=0.9,
    )
    for group in report["subtype"].values():
        assert group["observable_coverage_lower_ci"] == [0.5, 0.5]
        assert group["bootstrap_sampling"] == "within_primary_event_strata"


def test_paired_bootstrap_uses_the_same_stratified_selections():
    y = labels([10, 20, 30, 40])
    clinical = PredictionIntervals(np.zeros(4), np.array([15, 20, 25, 40]), 0.1, 100, 100)
    molecular = PredictionIntervals(np.zeros(4), clinical.upper + 10, 0.1, 100, 100)
    report = paired_interval_comparison(
        y,
        clinical,
        molecular,
        labels([150]),
        bootstrap_replicates=40,
        seed=3,
        bootstrap_strata=np.array([True, False, True, False]),
        confidence_level=0.95,
        min_censoring_survival=0.05,
        minimum_bootstrap_success_fraction=0.9,
    )
    assert report["bootstrap_sampling"] == "within_primary_event_strata"
    assert report["differences"]["width_median_months"]["ci"] == [10, 10]
    assert report["differences"]["width_median_months"]["status"] == "ok"


def test_bootstrap_rejects_misaligned_strata_and_invalid_success_gate():
    y = labels([10, 20])
    intervals = PredictionIntervals(np.zeros(2), np.full(2, 50), 0.1, 100, 100)
    with pytest.raises(ValueError, match="align"):
        coverage_summary(
            y,
            intervals,
            labels([150]),
            bootstrap_replicates=20,
            seed=1,
            bootstrap_strata=np.array([True]),
        )
    with pytest.raises(ValueError, match="minimum success fraction"):
        coverage_summary(
            y,
            intervals,
            labels([150]),
            bootstrap_replicates=20,
            seed=1,
            minimum_bootstrap_success_fraction=0,
        )
