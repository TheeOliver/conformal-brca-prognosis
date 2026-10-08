"""Synthetic contracts for one-sided restricted-time score-envelope calibration."""

import json

import numpy as np
import pytest

from brca.conformal.lower import (
    LowerBoundCalibration,
    calibrate_lower_bound,
    lower_score_bounds,
    signed_finite_sample_quantile,
)


def labels(time, event=None):
    time = np.asarray(time, dtype=float)
    result = np.empty(len(time), dtype=[("event", bool), ("time", float)])
    result["time"] = time
    result["event"] = True if event is None else event
    return result


@pytest.mark.parametrize("alpha,expected", [(0.2, -20), (0.1, -10), (0.05, -5)])
def test_exact_signed_order_statistics(alpha, expected):
    scores = np.arange(-99, 0, dtype=float)
    assert signed_finite_sample_quantile(scores, alpha) == expected
    assert signed_finite_sample_quantile(scores[::-1], alpha) == expected


def test_signed_correction_raises_conservative_base_and_round_trips():
    state = calibrate_lower_bound(labels(np.arange(11, 20)), np.full(9, 10), 0.2, horizon=30)
    assert state.quantile == -2
    predicted = state.predict(np.array([10, 29, np.inf]))
    np.testing.assert_array_equal(predicted.lower, [12, 30, 30])
    np.testing.assert_array_equal(predicted.upper, [30, 30, 30])
    np.testing.assert_array_equal(predicted.width, [18, 0, 0])
    assert predicted.n_cal == 9
    assert predicted.method == "conservative_lower_bound"
    restored = LowerBoundCalibration.from_dict(
        json.loads(json.dumps(state.to_dict(), allow_nan=False))
    )
    np.testing.assert_array_equal(restored.predict(np.array([10])).lower, [12])


def test_finite_rank_empty_missing_and_undersized_strata_return_full_support():
    assert np.isinf(signed_finite_sample_quantile(np.array([]), 0.1))
    assert np.isinf(signed_finite_sample_quantile(np.arange(-9, 0), 0.05))
    assert signed_finite_sample_quantile(np.arange(-19, 0), 0.05) == -1
    empty = calibrate_lower_bound(labels([]), np.array([]), 0.1, horizon=100)
    assert empty.n_cal == 0
    assert empty.predict(np.array([np.inf])).width[0] == 100
    state = calibrate_lower_bound(
        labels(np.arange(1, 51)),
        np.full(50, 10),
        0.1,
        horizon=100,
        strata=np.ones(50, dtype=bool),
        expected_strata=[False, True],
    )
    assert state.stratum_counts == {"False": 0, "True": 50}
    assert np.isinf(state.quantile)
    restored = LowerBoundCalibration.from_dict(
        json.loads(json.dumps(state.to_dict(), allow_nan=False))
    )
    predicted = restored.predict(np.array([0, 20, np.inf]))
    np.testing.assert_array_equal(predicted.lower, [0, 0, 0])
    np.testing.assert_array_equal(predicted.width, [100, 100, 100])
    assert predicted.quantile_is_infinite


def test_worst_signed_stratum_threshold_is_used_without_test_strata():
    y = labels(np.r_[np.arange(11, 20), np.arange(21, 30)])
    state = calibrate_lower_bound(
        y,
        np.full(18, 10),
        0.2,
        horizon=50,
        strata=np.repeat([False, True], 9),
        expected_strata=[False, True],
    )
    assert state.stratum_quantiles == {"False": -2, "True": -12}
    assert state.quantile == -2
    assert state.n_cal == 18
    assert state.predict(np.array([10])).lower[0] == 12


def test_zeros_ties_and_unreached_quantiles_are_handled_before_subtraction():
    assert signed_finite_sample_quantile(np.zeros(19), 0.05) == 0
    y = labels([20, 20, 100, 150], [True, False, False, True])
    scores = lower_score_bounds(y, np.array([30, 30, np.inf, np.inf]), horizon=100)
    np.testing.assert_array_equal(scores, [10, 10, 0, 0])
    state = calibrate_lower_bound(labels(np.full(19, 20)), np.full(19, 20), 0.05, horizon=100)
    assert state.quantile == 0
    assert state.predict(np.array([20])).lower[0] == 20
    clipped = LowerBoundCalibration(alpha=0.1, n_cal=99, horizon=100, quantile=80)
    np.testing.assert_array_equal(clipped.predict(np.array([0, np.inf])).lower, [0, 20])


def test_observed_scores_dominate_latent_scores_under_arbitrary_censoring():
    latent = np.array([20, 50, 150, 100, 70, 10], dtype=float)
    censor = np.array([30, 10, 120, 100, 60, 1], dtype=float)
    base = np.array([30, 20, np.inf, 110, 0, 40], dtype=float)
    observed = labels(np.minimum(latent, censor), latent <= censor)
    envelope = lower_score_bounds(observed, base, horizon=100)
    truth = np.minimum(base, 100) - np.minimum(latent, 100)
    assert np.all(envelope >= truth)
    known = observed["event"] | (observed["time"] >= 100)
    np.testing.assert_array_equal(envelope[known], truth[known])


@pytest.mark.parametrize("alpha", [0.2, 0.1, 0.05])
def test_repeated_uncensored_populations_achieve_exact_rank_coverage(alpha):
    rng = np.random.default_rng(48317)
    n_cal, n_test, repeats = 199, 250, 200
    rates = []
    for _ in range(repeats):
        x = rng.uniform(0, 1, n_cal + n_test)
        latent = 20 + 30 * x + rng.exponential(20, len(x))
        base = 20 + 30 * x - 20 * np.log1p(-alpha)
        state = calibrate_lower_bound(labels(latent[:n_cal]), base[:n_cal], alpha, horizon=200)
        intervals = state.predict(base[n_cal:])
        target = np.minimum(latent[n_cal:], 200)
        rates.append(np.mean((intervals.lower <= target) & (target <= intervals.upper)))
    expected = np.ceil((n_cal + 1) * (1 - alpha)) / (n_cal + 1)
    # Calibration/test datasets, rather than correlated patients, are the MC units.
    mc_se = np.std(rates, ddof=1) / np.sqrt(repeats)
    assert abs(np.mean(rates) - expected) < 4 * mc_se + 0.001


@pytest.mark.parametrize("mode", ["independent", "covariate_dependent", "informative"])
@pytest.mark.parametrize("alpha", [0.2, 0.1, 0.05])
def test_repeated_censored_populations_keep_conservative_coverage(mode, alpha):
    rng = np.random.default_rng(83115)
    n_cal, n_test, repeats = 199, 250, 150
    rates, nonzero_lower = [], []
    for _ in range(repeats):
        x = rng.uniform(0, 1, n_cal + n_test)
        latent = 30 + 25 * x + rng.exponential(20, len(x))
        base = 30 + 25 * x - 20 * np.log1p(-alpha)
        if mode == "independent":
            censor = rng.uniform(40, 130, len(x))
        elif mode == "covariate_dependent":
            censor = 45 + 25 * x + rng.uniform(5, 60, len(x))
        else:
            censor = latent * rng.uniform(0.75, 1.25, len(x))
        y = labels(np.minimum(latent, censor), latent <= censor)
        state = calibrate_lower_bound(y[:n_cal], base[:n_cal], alpha, horizon=120)
        intervals = state.predict(base[n_cal:])
        target = np.minimum(latent[n_cal:], 120)
        rates.append(np.mean((intervals.lower <= target) & (target <= intervals.upper)))
        nonzero_lower.append(np.mean(intervals.lower > 0))
    mc_se = np.std(rates, ddof=1) / np.sqrt(repeats)
    assert np.mean(rates) >= 1 - alpha - 4 * mc_se - 0.001
    # Guard against a vacuous implementation returning [0,tau] for every patient.
    assert np.mean(nonzero_lower) > 0.9


def test_primary_event_stratified_informative_censoring_covers_each_stratum():
    rng = np.random.default_rng(14739)
    group_rates = {False: [], True: []}
    for _ in range(150):
        x = rng.uniform(0, 1, 2400)
        latent = 30 + 25 * x + rng.exponential(20, len(x))
        censor = latent * rng.uniform(0.75, 1.25, len(x))
        y = labels(np.minimum(latent, censor), latent <= censor)
        base = 30 + 25 * x - 20 * np.log1p(-0.1)
        cal_indices, test_indices = [], []
        for group in [False, True]:
            indices = np.flatnonzero(y["event"] == group)
            rng.shuffle(indices)
            cal_indices.extend(indices[:99])
            test_indices.extend(indices[99:199])
        cal_indices, test_indices = np.array(cal_indices), np.array(test_indices)
        state = calibrate_lower_bound(
            y[cal_indices],
            base[cal_indices],
            0.1,
            horizon=120,
            strata=y["event"][cal_indices],
            expected_strata=[False, True],
        )
        intervals = state.predict(base[test_indices])
        target = np.minimum(latent[test_indices], 120)
        covered = (intervals.lower <= target) & (target <= intervals.upper)
        for group in [False, True]:
            group_rates[group].append(np.mean(covered[y["event"][test_indices] == group]))
    for rates in group_rates.values():
        mc_se = np.std(rates, ddof=1) / np.sqrt(len(rates))
        assert np.mean(rates) >= 0.9 - 4 * mc_se - 0.001


@pytest.mark.parametrize("horizon", [0, -1, np.nan, np.inf, -np.inf])
def test_nonfinite_or_nonpositive_horizon_is_rejected(horizon):
    with pytest.raises(ValueError, match="horizon"):
        calibrate_lower_bound(labels([10]), np.array([2]), 0.1, horizon=horizon)


@pytest.mark.parametrize("alpha", [0, 1, -0.1, np.nan])
def test_invalid_alpha_is_rejected(alpha):
    with pytest.raises(ValueError, match="alpha"):
        signed_finite_sample_quantile(np.array([-1, 0]), alpha)


@pytest.mark.parametrize("base", [[np.nan], [-1], [-np.inf], [[1]]])
def test_invalid_base_is_rejected(base):
    with pytest.raises(ValueError, match="Base lower"):
        lower_score_bounds(labels([10]), np.array(base), horizon=100)


def test_misaligned_inputs_and_unknown_strata_are_rejected():
    with pytest.raises(ValueError, match="aligned"):
        calibrate_lower_bound(labels([10, 20]), np.array([2]), 0.1, horizon=100)
    with pytest.raises(ValueError, match="expected_strata"):
        calibrate_lower_bound(
            labels([10]), np.array([2]), 0.1, horizon=100, expected_strata=[False, True]
        )
    with pytest.raises(ValueError, match="expected_strata"):
        calibrate_lower_bound(
            labels([10]),
            np.array([2]),
            0.1,
            horizon=100,
            strata=["unknown"],
            expected_strata=[False, True],
        )
