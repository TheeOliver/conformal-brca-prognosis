"""Numerically stabilized, explicitly exploratory Qin Cox bootstrap sensitivity.

The original curve and every bootstrap fit use the same fixed positive ridge
penalty. Every requested resample must succeed: no success-conditioned pivot
distribution is returned. Marginal censoring assumptions and weak event-pool
support remain limitations; ridge supplies no new coverage guarantee.
"""

from __future__ import annotations

import warnings
from collections import Counter
from numbers import Real

import numpy as np
import pandas as pd
from scipy.linalg import LinAlgWarning
from sklearn.exceptions import ConvergenceWarning
from sksurv.linear_model import CoxPHSurvivalAnalysis

from brca.conformal.qin import QinBootstrap, event_pool_weights
from brca.models.base import (
    config_value,
    step_survival_on_grid,
    validate_features,
    validate_survival_labels,
)
from brca.models.cox import CoxPHModel


class QinRidgeBootstrapError(ValueError):
    """The complete planned bootstrap failed its strict success requirement."""

    def __init__(self, diagnostics: dict):
        self.diagnostics = diagnostics
        super().__init__(
            "Qin ridge requires all bootstrap refits to succeed: "
            f"{diagnostics['bootstrap_successful']}/{diagnostics['bootstrap_requested']}"
        )


def _ridge_parameters(cox_config) -> dict:
    parameters = {
        name: config_value(cox_config, name) for name in ("alpha", "ties", "n_iter", "tol")
    }
    alpha = parameters["alpha"]
    if (
        isinstance(alpha, bool)
        or not isinstance(alpha, Real)
        or not np.isfinite(alpha)
        or alpha <= 0
    ):
        raise ValueError("Qin ridge requires a fixed positive finite scalar alpha")
    if parameters["ties"] not in {"breslow", "efron"}:
        raise ValueError("Cox ties must be breslow or efron")
    n_iter, tol = parameters["n_iter"], parameters["tol"]
    if isinstance(n_iter, bool) or not isinstance(n_iter, int) or n_iter < 1:
        raise ValueError("Cox n_iter must be a positive integer")
    if not isinstance(tol, Real) or not np.isfinite(tol) or tol <= 0:
        raise ValueError("Cox tol must be positive and finite")
    return parameters


def _require_event_time_fields(y: np.ndarray, n: int) -> None:
    validate_survival_labels(y, n)
    if y.dtype.names != ("event", "time"):
        raise ValueError("Qin requires structured event/time survival fields")


def _raise_numerical_warnings() -> None:
    for warning in (ConvergenceWarning, LinAlgWarning, RuntimeWarning):
        warnings.simplefilter("error", warning)


def fit_qin_ridge_model(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    cox_config,
    *,
    rng: np.random.Generator,
) -> CoxPHModel:
    """Fit the original ridge curve on an already prepared training matrix.

    Uses the common Cox adapter, which omits original-training constants and
    reports that its unpenalized PH score test is unavailable for ridge fits.
    Training preprocessing stays fixed across subsequent bootstrap refits.
    """
    _ridge_parameters(cox_config)
    validate_features(X_train)
    _require_event_time_fields(y_train, len(X_train))
    with warnings.catch_warnings():
        _raise_numerical_warnings()
        model = CoxPHModel(cox_config, rng=rng).fit(X_train, y_train)
    return model


def bootstrap_ridge_cox_pivots(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    original_model: CoxPHModel,
    cox_config,
    *,
    replicates: int,
    minimum_censoring_survival: float,
    rng: np.random.Generator,
) -> QinBootstrap:
    """Run all planned ridge refits; any failure rejects the complete result.

    Sampling order matches the original Qin implementation: draw the training
    bootstrap indices, then one inverse-censoring-weighted original event.
    Constant and rank-deficient bootstrap matrices are retained. Ridge gives
    their coefficients a defined shrinkage solution, including predictions for
    original event subjects whose category disappears from a resample.
    """
    parameters = _ridge_parameters(cox_config)
    validate_features(X_train, original_model.feature_names_)
    _require_event_time_fields(y_train, len(X_train))
    if isinstance(replicates, bool) or not isinstance(replicates, int) or replicates < 2:
        raise ValueError("At least two integer bootstrap replicates are required")
    original_parameters = original_model.estimator_.get_params(deep=False)
    for name, value in parameters.items():
        if not np.array_equal(np.asarray(original_parameters[name]), np.asarray(value)):
            raise ValueError(f"Original and bootstrap Cox parameters differ: {name}")
    active = list(original_model.active_feature_names_)
    if not active or not set(active).issubset(X_train.columns):
        raise ValueError("Original Cox active feature schema is invalid")
    if tuple(original_model.estimator_.feature_names_in_) != tuple(active):
        raise ValueError("Original Cox estimator and active feature schema differ")
    if not np.isfinite(original_model.estimator_.coef_).all():
        raise ValueError("Original Cox coefficients are nonfinite")

    event_indices, probabilities, event_diagnostics = event_pool_weights(
        y_train, minimum_censoring_survival=minimum_censoring_survival
    )
    X = X_train.loc[:, active]
    values = X.to_numpy(float)
    pivots = []
    failures = Counter()
    constant_counts = Counter()
    constant_replicates = rank_deficient_replicates = 0
    failed_constant_replicates = failed_rank_deficient_replicates = 0
    no_event_replicates = tail_carry_replicates = 0
    for _ in range(replicates):
        sampled = rng.integers(0, len(X), size=len(X))
        event_index = int(rng.choice(event_indices, p=probabilities))
        sampled_values, sampled_y = values[sampled], y_train[sampled]
        constant = np.ptp(sampled_values, axis=0) == 0
        rank_deficient = np.linalg.matrix_rank(sampled_values - sampled_values.mean(axis=0)) < len(
            active
        )
        has_constant = bool(constant.any())
        constant_replicates += int(has_constant)
        rank_deficient_replicates += int(rank_deficient)
        constant_counts.update(np.asarray(active)[constant].tolist())
        has_event = bool(sampled_y["event"].any())
        no_event_replicates += int(not has_event)
        failure_reason = None
        try:
            if not has_event:
                failure_reason = "no_observed_events"
                raise ValueError("No observed event in bootstrap sample")
            estimator = CoxPHSurvivalAnalysis(**parameters)
            with warnings.catch_warnings():
                _raise_numerical_warnings()
                estimator.fit(X.iloc[sampled], sampled_y)
                if not np.isfinite(estimator.coef_).all():
                    failure_reason = "nonfinite_coefficients"
                    raise ValueError("Nonfinite bootstrap coefficients")
                survival = estimator.predict_survival_function(
                    X.iloc[[event_index]], return_array=True
                )
                event_time = float(y_train["time"][event_index])
                pivot = float(
                    step_survival_on_grid(estimator.unique_times_, survival, [event_time])[0, 0]
                )
                if not np.isfinite(pivot) or not 0 <= pivot <= 1:
                    failure_reason = "invalid_survival_pivot"
                    raise ValueError("Invalid bootstrap survival pivot")
            pivots.append(pivot)
            tail_carry_replicates += int(event_time > np.max(sampled_y["time"][sampled_y["event"]]))
        except (
            ArithmeticError,
            ConvergenceWarning,
            LinAlgWarning,
            ValueError,
            RuntimeError,
            RuntimeWarning,
            np.linalg.LinAlgError,
        ) as error:
            failures[failure_reason or type(error).__name__] += 1
            failed_constant_replicates += int(has_constant)
            failed_rank_deficient_replicates += int(rank_deficient)

    diagnostics = {
        **event_diagnostics,
        "bootstrap_requested": replicates,
        "bootstrap_attempted": replicates,
        "bootstrap_successful": len(pivots),
        "bootstrap_failed": sum(failures.values()),
        "bootstrap_success_fraction": len(pivots) / replicates,
        "required_bootstrap_success_fraction": 1.0,
        "failure_reasons": dict(sorted(failures.items())),
        "bootstrap_tail_carry_replicates": tail_carry_replicates,
        "bootstrap_design": {
            "active_feature_count": len(active),
            "replicates_with_constant_columns": constant_replicates,
            "replicates_with_rank_deficiency": rank_deficient_replicates,
            "constant_column_replicate_counts": dict(sorted(constant_counts.items())),
            "failed_replicates_with_constant_columns": failed_constant_replicates,
            "failed_replicates_with_rank_deficiency": failed_rank_deficient_replicates,
            "replicates_without_events": no_event_replicates,
        },
        "working_model": "fixed_ridge_cox_ph",
        "working_model_parameters": parameters,
        "original_bootstrap_parameter_match": True,
        "preprocessing": "original_training_fit_fixed_across_bootstrap",
        "censoring_assumption": "marginal independent censoring C independent of (T,X)",
        "coverage_guarantee": "exploratory_ridge_adaptation_no_new_guarantee",
    }
    if len(pivots) != replicates:
        raise QinRidgeBootstrapError(diagnostics)
    return QinBootstrap(np.asarray(pivots), diagnostics)
