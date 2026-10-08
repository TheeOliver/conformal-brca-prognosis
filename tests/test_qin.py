"""Synthetic contract tests for the exploratory Qin Cox comparator."""

import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from sklearn.exceptions import ConvergenceWarning
from sksurv.util import Surv

from brca.conformal.qin import (
    QinBootstrap,
    bootstrap_cox_pivots,
    event_pool_weights,
    inverse_survival_times,
    predict_qin_intervals,
)


def test_qin_event_pool_uses_reverse_km_inverse_weights():
    y = Surv.from_arrays([True, False, True, False], [1.0, 2.0, 3.0, 4.0])
    indices, probabilities, diagnostics = event_pool_weights(y, minimum_censoring_survival=0.0)
    np.testing.assert_array_equal(indices, [0, 2])
    np.testing.assert_allclose(probabilities, [0.4, 0.6])
    assert diagnostics["n_training_events"] == 2
    with pytest.raises(ValueError, match="censoring support"):
        event_pool_weights(y, minimum_censoring_survival=0.7)


def test_qin_event_before_censor_tie_convention():
    y = Surv.from_arrays([True, False, True], [1.0, 1.0, 2.0])
    indices, probabilities, _ = event_pool_weights(y, minimum_censoring_survival=0.0)
    np.testing.assert_array_equal(indices, [0, 2])
    np.testing.assert_allclose(probabilities, [1 / 3, 2 / 3])


class _ExponentialStepEstimator:
    def __init__(self):
        self.unique_times_ = np.linspace(0.01, 20.0, 2000)

    def predict_survival_function(self, X, *, return_array):
        assert return_array
        return np.exp(-self.unique_times_[None, :]) * np.ones((len(X), 1))


def test_qin_pivot_inversion_empirical_coverage_on_uncensored_synthetic_data():
    rng = np.random.default_rng(7)
    bootstrap = QinBootstrap(rng.uniform(size=5000), {"synthetic": True})
    model = SimpleNamespace(estimator_=_ExponentialStepEstimator(), active_feature_names_=("x",))
    X_test = pd.DataFrame({"x": np.zeros(10000)})
    intervals = predict_qin_intervals(model, X_test, bootstrap, alpha=0.1, horizon=20.0)
    true_times = rng.exponential(size=len(X_test))
    coverage = np.mean((true_times >= intervals.lower) & (true_times <= intervals.upper))
    assert 0.87 <= coverage <= 0.93
    assert 0 < np.median(intervals.width) < intervals.horizon
    assert intervals.n_cal == 5000


def test_qin_upper_endpoint_includes_survival_plateau_at_lower_pivot():
    survival = np.array([[0.9, 0.5, 0.5, 0.2]])
    times = np.array([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(inverse_survival_times(survival, times, 0.5), [2.0])
    np.testing.assert_array_equal(inverse_survival_times(survival, times, 0.5, strict=True), [4.0])
    np.testing.assert_array_equal(inverse_survival_times(survival, times, 1.0), [0.0])
    np.testing.assert_array_equal(inverse_survival_times(survival, times, 1.0, strict=True), [1.0])


def test_qin_bootstrap_refits_synthetic_cox_and_counts_nonconvergence(monkeypatch):
    rng = np.random.default_rng(81)
    X = pd.DataFrame({"x": rng.normal(size=80)})
    event_times = rng.exponential(scale=8, size=80)
    censor_times = rng.exponential(scale=20, size=80)
    y = Surv.from_arrays(event_times <= censor_times, np.minimum(event_times, censor_times))
    from sksurv.linear_model import CoxPHSurvivalAnalysis

    estimator = CoxPHSurvivalAnalysis(alpha=0.0, ties="breslow", n_iter=100, tol=1e-9).fit(X, y)
    model = SimpleNamespace(
        estimator_=estimator,
        active_feature_names_=("x",),
        feature_names_=("x",),
    )
    cfg = SimpleNamespace(alpha=0.0, ties="breslow", n_iter=100, tol=1e-9)
    original_fit = CoxPHSurvivalAnalysis.fit
    calls = 0

    def warn_once(self, features, labels):
        nonlocal calls
        calls += 1
        if calls == 1:
            warnings.warn("synthetic nonconvergence", ConvergenceWarning, stacklevel=2)
        return original_fit(self, features, labels)

    monkeypatch.setattr(CoxPHSurvivalAnalysis, "fit", warn_once)
    result = bootstrap_cox_pivots(
        X,
        y,
        model,
        cfg,
        replicates=30,
        minimum_success_fraction=0.9,
        minimum_censoring_survival=0.0,
        rng=np.random.default_rng(92),
    )
    assert len(result.pivots) == 29
    assert np.all((result.pivots >= 0) & (result.pivots <= 1))
    assert result.diagnostics["bootstrap_requested"] == 30
    assert result.diagnostics["bootstrap_failed"] == 1


@pytest.mark.slow
def test_qin_end_to_end_coverage_with_independent_censoring():
    """A seeded synthetic survival population exercises weighting, refits and inversion."""
    from sksurv.linear_model import CoxPHSurvivalAnalysis

    rng = np.random.default_rng(481)
    n_train, n_test = 500, 4000
    x_train, x_test = rng.normal(size=n_train), rng.normal(size=n_test)
    train_event = rng.exponential(scale=10 / np.exp(0.5 * x_train))
    train_censor = rng.exponential(scale=35, size=n_train)
    y_train = Surv.from_arrays(train_event <= train_censor, np.minimum(train_event, train_censor))
    X_train = pd.DataFrame({"x": x_train})
    X_test = pd.DataFrame({"x": x_test})
    estimator = CoxPHSurvivalAnalysis(alpha=0.0, ties="breslow", n_iter=100, tol=1e-9).fit(
        X_train, y_train
    )
    model = SimpleNamespace(
        estimator_=estimator,
        active_feature_names_=("x",),
        feature_names_=("x",),
    )
    config = SimpleNamespace(alpha=0.0, ties="breslow", n_iter=100, tol=1e-9)
    bootstrap = bootstrap_cox_pivots(
        X_train,
        y_train,
        model,
        config,
        replicates=400,
        minimum_success_fraction=0.9,
        minimum_censoring_survival=0.0,
        rng=np.random.default_rng(482),
    )
    intervals = predict_qin_intervals(model, X_test, bootstrap, alpha=0.1, horizon=80.0)
    test_event = rng.exponential(scale=10 / np.exp(0.5 * x_test))
    target = np.minimum(test_event, 80.0)
    coverage = np.mean((target >= intervals.lower) & (target <= intervals.upper))
    assert 0.85 <= coverage <= 0.95
    assert 0 < np.median(intervals.width) < 80
