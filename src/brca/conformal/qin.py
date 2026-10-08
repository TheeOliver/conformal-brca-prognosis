"""Qin et al. (2025) bootstrap CPI, Cox working-model implementation.

Uses marginal reverse-KM censoring weights on observed training events and
bootstrap refits of a Cox model. This is an approximate resampling CPI, not a
split-conformal finite-sample guarantee. Preprocessing is frozen from train.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.linalg import LinAlgWarning
from sklearn.exceptions import ConvergenceWarning
from sksurv.linear_model import CoxPHSurvivalAnalysis

from brca.conformal.survival import CensoringDistribution, PredictionIntervals
from brca.models.base import step_survival_on_grid, validate_survival_labels


@dataclass(frozen=True)
class QinBootstrap:
    pivots: np.ndarray
    diagnostics: dict

    def thresholds(self, alpha: float) -> tuple[float, float]:
        if not 0 < alpha < 1:
            raise ValueError("alpha must be miscoverage in (0, 1)")
        if len(self.pivots) < 2:
            raise ValueError("At least two bootstrap pivots are required")
        lower, upper = np.quantile(self.pivots, [alpha / 2, 1 - alpha / 2])
        return float(lower), float(upper)


def event_pool_weights(
    y_train: np.ndarray, *, minimum_censoring_survival: float
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Return event indices and probabilities proportional to 1/G(Y-)."""
    validate_survival_labels(y_train, len(y_train))
    if not 0 <= minimum_censoring_survival < 1:
        raise ValueError("Invalid censoring support threshold")
    indices = np.flatnonzero(y_train["event"])
    censoring = CensoringDistribution.fit(y_train)
    support = censoring.survival_ge(y_train["time"][indices])
    if not np.isfinite(support).all() or np.any(support <= minimum_censoring_survival):
        raise ValueError("Training event times lack the required censoring support")
    raw = 1 / support
    probabilities = raw / raw.sum()
    return (
        indices,
        probabilities,
        {
            "n_training_events": len(indices),
            "minimum_event_censoring_survival": float(support.min()),
            "maximum_event_weight": float(raw.max()),
            "event_pool_effective_sample_size": float(raw.sum() ** 2 / np.square(raw).sum()),
        },
    )


def bootstrap_cox_pivots(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    original_model,
    cox_config,
    *,
    replicates: int,
    minimum_success_fraction: float,
    minimum_censoring_survival: float,
    rng: np.random.Generator,
) -> QinBootstrap:
    """Bootstrap Cox refit, then draw one IPCW event from the original train.

    The saved original model and train-fitted preprocessing are fixed. Each
    bootstrap sample refits Cox coefficients and baseline survival. Nonconverged
    replicates are discarded only if the prespecified success gate passes.
    """
    validate_survival_labels(y_train, len(X_train))
    if not isinstance(replicates, int) or replicates < 2 or not 0 < minimum_success_fraction <= 1:
        raise ValueError("Invalid bootstrap count or success fraction")
    if tuple(X_train.columns) != original_model.feature_names_:
        raise ValueError("Training features differ from the saved Cox model")
    event_indices, probabilities, diagnostics = event_pool_weights(
        y_train, minimum_censoring_survival=minimum_censoring_survival
    )
    active = list(original_model.active_feature_names_)
    X = X_train.loc[:, active]
    pivots = []
    failures = 0
    tail_extrapolations = 0
    for _ in range(replicates):
        sampled = rng.integers(0, len(X), size=len(X))
        event_index = int(rng.choice(event_indices, p=probabilities))
        estimator = CoxPHSurvivalAnalysis(
            alpha=cox_config.alpha,
            ties=cox_config.ties,
            n_iter=cox_config.n_iter,
            tol=cox_config.tol,
        )
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", ConvergenceWarning)
                warnings.simplefilter("error", LinAlgWarning)
                warnings.simplefilter("error", RuntimeWarning)
                estimator.fit(X.iloc[sampled], y_train[sampled])
            if not np.isfinite(estimator.coef_).all():
                raise ValueError("Nonfinite bootstrap coefficients")
            survival = estimator.predict_survival_function(X.iloc[[event_index]], return_array=True)
            event_time = float(y_train["time"][event_index])
            if event_time > np.max(y_train["time"][sampled][y_train["event"][sampled]]):
                tail_extrapolations += 1
            pivot = float(
                step_survival_on_grid(estimator.unique_times_, survival, [event_time])[0, 0]
            )
            if not np.isfinite(pivot) or not 0 <= pivot <= 1:
                raise ValueError("Invalid bootstrap survival pivot")
            pivots.append(pivot)
        except (
            ArithmeticError,
            ConvergenceWarning,
            LinAlgWarning,
            ValueError,
            RuntimeError,
            RuntimeWarning,
            np.linalg.LinAlgError,
        ):
            failures += 1
    if len(pivots) / replicates < minimum_success_fraction or len(pivots) < 2:
        raise ValueError(f"Qin bootstrap success gate failed: {len(pivots)}/{replicates}")
    return QinBootstrap(
        np.asarray(pivots),
        {
            **diagnostics,
            "bootstrap_requested": replicates,
            "bootstrap_successful": len(pivots),
            "bootstrap_failed": failures,
            "bootstrap_tail_carry_replicates": tail_extrapolations,
            "bootstrap_success_fraction": len(pivots) / replicates,
            "working_model": "cox_ph",
            "preprocessing": "original_training_fit_fixed_across_bootstrap",
            "censoring_assumption": "marginal independent censoring C independent of (T,X)",
            "coverage_guarantee": "approximate_asymptotic_under_assumptions",
        },
    )


def inverse_survival_times(
    survival: np.ndarray, times: np.ndarray, threshold: float, *, strict: bool = False
) -> np.ndarray:
    """First event knot where S<=threshold (or S<threshold for the upper edge)."""
    survival, times = np.asarray(survival, float), np.asarray(times, float)
    if survival.ndim != 2 or survival.shape[1] != len(times) or not 0 <= threshold <= 1:
        raise ValueError("Invalid survival grid or probability threshold")
    if np.any(np.diff(survival, axis=1) > 1e-10):
        raise ValueError("Survival curves must be nonincreasing")
    if threshold == 1 and not strict:
        return np.zeros(len(survival))
    crossed = survival < threshold if strict else survival <= threshold
    result = np.full(len(survival), np.inf)
    found = crossed.any(axis=1)
    result[found] = times[crossed[found].argmax(axis=1)]
    return result


def predict_qin_intervals(
    original_model, X_test: pd.DataFrame, bootstrap: QinBootstrap, *, alpha: float, horizon: float
) -> PredictionIntervals:
    """Invert original Cox S at upper/lower bootstrap pivot quantiles."""
    if not np.isfinite(horizon) or horizon <= 0:
        raise ValueError("Restricted target requires a positive finite horizon")
    lower_pivot, upper_pivot = bootstrap.thresholds(alpha)
    survival = original_model.estimator_.predict_survival_function(
        X_test.loc[:, original_model.active_feature_names_], return_array=True
    )
    knots = original_model.estimator_.unique_times_
    lower = inverse_survival_times(survival, knots, upper_pivot)
    upper = inverse_survival_times(survival, knots, lower_pivot, strict=True)
    return PredictionIntervals(
        np.minimum(lower, horizon),
        np.minimum(upper, horizon),
        alpha,
        len(bootstrap.pivots),
        horizon,
        "qin_bootstrap_cox",
    )
