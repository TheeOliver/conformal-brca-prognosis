"""Synthetic tests for strict fixed-ridge Qin bootstrap stabilization."""

import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from sklearn.exceptions import ConvergenceWarning
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv

from brca.conformal.qin import event_pool_weights, predict_qin_intervals
from brca.conformal.qin_ridge import (
    QinRidgeBootstrapError,
    bootstrap_ridge_cox_pivots,
    fit_qin_ridge_model,
)


def _config(**updates):
    return SimpleNamespace(
        **{"alpha": 1.0, "ties": "breslow", "n_iter": 100, "tol": 1e-9, **updates}
    )


def _rare_indicator_training():
    rng = np.random.default_rng(711)
    n = 70
    X = pd.DataFrame({"x": rng.normal(size=n), "rare": np.arange(n) == 0}).astype(float)
    times = rng.exponential(scale=8 / np.exp(0.4 * X["x"].to_numpy()))
    censoring = rng.exponential(scale=30, size=n)
    # The sole rare-category subject is observed, and can be drawn as a pivot
    # even when its category disappears from the independently drawn bootstrap.
    censoring[0] = times[0] + 1
    y = Surv.from_arrays(times <= censoring, np.minimum(times, censoring))
    return X, y


def _fit(X, y, config=None):
    return fit_qin_ridge_model(X, y, config or _config(), rng=np.random.default_rng(712))


def _bootstrap(X, y, model, *, config=None, replicates=35, seed=713):
    return bootstrap_ridge_cox_pivots(
        X,
        y,
        model,
        config or _config(),
        replicates=replicates,
        minimum_censoring_survival=0.0,
        rng=np.random.default_rng(seed),
    )


def test_ridge_handles_absent_category_and_reports_exact_design_counts():
    X, y = _rare_indicator_training()
    model = _fit(X, y)
    bootstrap = _bootstrap(X, y, model)
    event_indices, probabilities, _ = event_pool_weights(y, minimum_censoring_survival=0.0)
    rng = np.random.default_rng(713)
    absent_count = 0
    for _ in range(35):
        sampled = rng.integers(0, len(X), size=len(X))
        rng.choice(event_indices, p=probabilities)
        absent_count += int(not np.any(sampled == 0))
    assert absent_count > 0
    diagnostics = bootstrap.diagnostics
    assert diagnostics["bootstrap_successful"] == diagnostics["bootstrap_attempted"] == 35
    assert diagnostics["bootstrap_failed"] == 0
    assert diagnostics["failure_reasons"] == {}
    assert diagnostics["bootstrap_design"]["replicates_with_constant_columns"] == absent_count
    assert diagnostics["bootstrap_design"]["replicates_with_rank_deficiency"] == absent_count
    assert diagnostics["bootstrap_design"]["constant_column_replicate_counts"] == {
        "rare": absent_count
    }
    assert np.isfinite(bootstrap.pivots).all()
    assert np.all((bootstrap.pivots >= 0) & (bootstrap.pivots <= 1))
    assert model.ph_diagnostic_["reason"] == "ph_score_inference_requires_unpenalized_fit"


def test_ridge_predicts_pivot_subject_whose_category_is_absent(monkeypatch):
    X, y = _rare_indicator_training()
    model = _fit(X, y)
    original_predict = CoxPHSurvivalAnalysis.predict_survival_function
    predictions = 0

    def check_absent_category(self, features, **kwargs):
        nonlocal predictions
        assert features["rare"].iloc[0] == 1
        # A zero column has an exact zero ridge coefficient; no column is dropped.
        assert self.coef_[list(features.columns).index("rare")] == 0
        predictions += 1
        return original_predict(self, features, **kwargs)

    class PrescribedDraws:
        def integers(self, low, high, *, size):
            assert (low, high, size) == (0, len(X), len(X))
            return np.resize(np.arange(1, len(X)), size)

        def choice(self, choices, *, p):
            assert 0 in choices and np.isclose(p.sum(), 1)
            return 0

    monkeypatch.setattr(CoxPHSurvivalAnalysis, "predict_survival_function", check_absent_category)
    bootstrap = bootstrap_ridge_cox_pivots(
        X,
        y,
        model,
        _config(),
        replicates=2,
        minimum_censoring_survival=0.0,
        rng=PrescribedDraws(),
    )
    assert predictions == 2
    assert len(bootstrap.pivots) == 2


def test_ridge_handles_separated_event_pattern_without_failed_bootstraps():
    rng = np.random.default_rng(714)
    group = np.repeat([0.0, 1.0], 40)
    X = pd.DataFrame({"x": rng.normal(size=80), "group": group})
    # All group-1 events precede every group-0 censoring time: monotone
    # unpenalized partial likelihood in the group coefficient.
    times = np.where(group == 1, rng.uniform(0.5, 2, size=80), rng.uniform(4, 8, size=80))
    y = Surv.from_arrays(group == 1, times)
    model = _fit(X, y)
    bootstrap = _bootstrap(X, y, model, replicates=25)
    assert np.isfinite(model.estimator_.coef_).all()
    assert bootstrap.diagnostics["bootstrap_success_fraction"] == 1
    assert bootstrap.diagnostics["bootstrap_failed"] == 0


@pytest.mark.parametrize(
    "changed", [{"alpha": 2.0}, {"ties": "efron"}, {"n_iter": 200}, {"tol": 1e-8}]
)
def test_ridge_rejects_original_bootstrap_parameter_mismatch(changed):
    X, y = _rare_indicator_training()
    model = _fit(X, y)
    with pytest.raises(ValueError, match="Original and bootstrap Cox parameters differ"):
        _bootstrap(X, y, model, config=_config(**changed), replicates=2)


@pytest.mark.parametrize("alpha", [0.0, -1.0, np.inf, np.nan, [1.0, 1.0], True])
def test_ridge_requires_fixed_positive_scalar_penalty(alpha):
    X, y = _rare_indicator_training()
    with pytest.raises(ValueError, match="fixed positive finite scalar alpha"):
        _fit(X, y, _config(alpha=alpha))


def test_ridge_rejects_failed_bootstrap_instead_of_conditioning_on_success(monkeypatch):
    X, y = _rare_indicator_training()
    model = _fit(X, y)
    original_fit = CoxPHSurvivalAnalysis.fit
    calls = 0

    def warn_once(self, features, labels):
        nonlocal calls
        calls += 1
        if calls == 1:
            warnings.warn("synthetic nonconvergence", ConvergenceWarning, stacklevel=2)
        return original_fit(self, features, labels)

    monkeypatch.setattr(CoxPHSurvivalAnalysis, "fit", warn_once)
    with pytest.raises(QinRidgeBootstrapError) as caught:
        _bootstrap(X, y, model, replicates=12)
    assert calls == 12
    diagnostics = caught.value.diagnostics
    assert diagnostics["bootstrap_requested"] == diagnostics["bootstrap_attempted"] == 12
    assert diagnostics["bootstrap_successful"] == 11
    assert diagnostics["bootstrap_failed"] == 1
    assert diagnostics["failure_reasons"] == {"ConvergenceWarning": 1}
    assert diagnostics["required_bootstrap_success_fraction"] == 1


def test_ridge_reproducible_with_same_draws_and_original_fit():
    X, y = _rare_indicator_training()
    model = _fit(X, y)
    first = _bootstrap(X, y, model, replicates=15)
    second = _bootstrap(X, y, model, replicates=15)
    np.testing.assert_array_equal(first.pivots, second.pivots)
    assert first.diagnostics == second.diagnostics


@pytest.mark.slow
def test_ridge_end_to_end_coverage_across_synthetic_censored_training_populations():
    """Independent training populations exercise refits, weights and inversion."""
    coverages, widths = [], []
    for seed in (720, 721, 722):
        rng = np.random.default_rng(seed)
        n_train, n_test = 500, 3000
        x_train, x_test = rng.normal(size=n_train), rng.normal(size=n_test)
        true_train = rng.exponential(scale=10 / np.exp(0.5 * x_train))
        censoring = rng.exponential(scale=45, size=n_train)
        y = Surv.from_arrays(true_train <= censoring, np.minimum(true_train, censoring))
        X = pd.DataFrame({"x": x_train})
        model = _fit(X, y)
        bootstrap = _bootstrap(X, y, model, replicates=400, seed=seed + 100)
        intervals = predict_qin_intervals(
            model,
            pd.DataFrame({"x": x_test}),
            bootstrap,
            alpha=0.1,
            horizon=80.0,
        )
        target = np.minimum(rng.exponential(scale=10 / np.exp(0.5 * x_test)), 80.0)
        coverages.append(np.mean((target >= intervals.lower) & (target <= intervals.upper)))
        widths.append(np.median(intervals.width))
    # A regression check for approximate marginal calibration, not evidence of
    # a finite-sample theorem: three populations retain noticeable Monte Carlo error.
    assert 0.85 <= np.mean(coverages) <= 0.95
    assert np.min(coverages) >= 0.8
    assert np.all((np.asarray(widths) > 0) & (np.asarray(widths) < 80))
