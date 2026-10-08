"""Cox proportional hazards adapter and a training-only PH score diagnostic."""

from __future__ import annotations

import warnings
from typing import Any, Self

import numpy as np
import pandas as pd
from scipy.stats import chi2, rankdata
from sklearn.exceptions import ConvergenceWarning
from sksurv.linear_model import CoxPHSurvivalAnalysis

from brca.models.base import (
    StepSurvivalModel,
    config_value,
    validate_features,
    validate_survival_labels,
)


def proportional_hazards_test(
    X: pd.DataFrame,
    y: np.ndarray,
    linear_predictor: np.ndarray,
    *,
    ties: str,
    transform: str,
    significance_level: float,
) -> dict[str, Any]:
    """Score test for beta(t)=beta+gamma*g(t), with beta treated as nuisance.

    Risk-set covariance blocks implement the time-interaction score test described
    in Therneau's survival vignette, section 'Computational details'. Both Breslow
    and Efron tie likelihoods are supported. Only unpenalized fits are inferentially
    supported; callers must enforce this. The chi-square reference is asymptotic.
    Per-column tests add one interaction at a time; encoded factors are not pooled.
    """
    validate_features(X)
    validate_survival_labels(y, len(X))
    if ties not in {"breslow", "efron"}:
        raise ValueError("PH diagnostic ties must be breslow or efron")
    if transform not in {"log", "identity", "rank"}:
        raise ValueError("PH diagnostic transform must be log, identity or rank")
    if not 0 < significance_level < 1:
        raise ValueError("PH diagnostic significance level must lie within (0, 1)")
    linear_predictor = np.asarray(linear_predictor, dtype=float)
    if linear_predictor.shape != (len(X),) or not np.isfinite(linear_predictor).all():
        raise ValueError("PH diagnostic requires finite training linear predictors aligned to X")

    event, time = (y[name] for name in y.dtype.names)
    event_times, event_counts = np.unique(time[event], return_counts=True)
    if transform == "log":
        transformed = np.log(event_times)
    elif transform == "identity":
        transformed = event_times.copy()
    else:
        repeated_ranks = rankdata(time[event])
        transformed = np.array([repeated_ranks[time[event] == t][0] for t in event_times])
    transformed -= np.average(transformed, weights=event_counts)

    values = X.to_numpy(dtype=float)
    p = values.shape[1]
    base_score, interaction_score = np.zeros(p), np.zeros(p)
    h_base, h_cross, h_interaction = (np.zeros((p, p)) for _ in range(3))
    for event_time, event_count, g in zip(event_times, event_counts, transformed, strict=True):
        at_risk = time >= event_time
        deaths = event & (time == event_time)
        x_risk = values[at_risk]
        weights = np.exp(linear_predictor[at_risk] - linear_predictor[at_risk].max())
        risk_zero = weights.sum()
        risk_one = weights @ x_risk
        risk_two = (x_risk.T * weights) @ x_risk
        # Efron removes a fraction of the tied failures at each substep.
        event_weights = weights[deaths[at_risk]]
        x_event = values[deaths]
        event_zero = event_weights.sum()
        event_one = event_weights @ x_event
        event_two = (x_event.T * event_weights) @ x_event
        residual = x_event.sum(axis=0)
        covariance = np.zeros((p, p))
        for k in range(event_count):
            fraction = k / event_count if ties == "efron" else 0
            denominator = risk_zero - fraction * event_zero
            mean = (risk_one - fraction * event_one) / denominator
            residual -= mean
            covariance += (risk_two - fraction * event_two) / denominator - np.outer(mean, mean)
        base_score += residual
        interaction_score += g * residual
        h_base += covariance
        h_cross += g * covariance
        h_interaction += g * g * covariance

    result: dict[str, Any] = {
        "method": "time_interaction_score_test",
        "reference_distribution": "asymptotic_chi_square",
        "training_only": True,
        "ties": ties,
        "time_transform": transform,
        "n_training": len(X),
        "n_events": int(event.sum()),
        "significance_level": float(significance_level),
        "multiplicity_adjusted": False,
        "factor_tests": "individual_encoded_columns",
    }
    if np.linalg.matrix_rank(h_base) < p:
        return {**result, "status": "unavailable", "reason": "singular_base_information"}
    projected_score = interaction_score - h_cross.T @ np.linalg.solve(h_base, base_score)
    information = h_interaction - h_cross.T @ np.linalg.solve(h_base, h_cross)
    information = (information + information.T) / 2
    if np.linalg.matrix_rank(information) < p or np.linalg.eigvalsh(information).min() <= 0:
        return {**result, "status": "unavailable", "reason": "singular_interaction_information"}
    statistic = float(projected_score @ np.linalg.solve(information, projected_score))
    global_p = float(chi2.sf(statistic, p))
    columns = {}
    for j, name in enumerate(X.columns):
        column_statistic = float(projected_score[j] ** 2 / information[j, j])
        p_value = float(chi2.sf(column_statistic, 1))
        columns[name] = {
            "statistic": column_statistic,
            "degrees_of_freedom": 1,
            "p_value": p_value,
            "flagged": p_value < significance_level,
        }
    return {
        **result,
        "status": "ok",
        "global": {
            "statistic": statistic,
            "degrees_of_freedom": p,
            "p_value": global_p,
            "flagged": global_p < significance_level,
        },
        "columns": columns,
    }


class CoxPHModel(StepSurvivalModel):
    """Train on a prepared matrix; store an aggregate PH diagnostic with the fit."""

    def __init__(self, config: Any, *, rng: np.random.Generator):
        self.config = config
        self.rng = rng  # Cox optimization itself is deterministic.

    def fit(self, X: pd.DataFrame, y: np.ndarray) -> Self:
        validate_features(X)
        validate_survival_labels(y, len(X))
        active_features = tuple(X.columns[X.nunique(dropna=False) > 1])
        if not active_features:
            raise ValueError("Cox requires at least one nonconstant training feature")
        active_X = X.loc[:, active_features]
        if np.linalg.matrix_rank(active_X.to_numpy() - active_X.to_numpy().mean(axis=0)) < len(
            active_features
        ):
            raise ValueError("Cox training features are collinear; omit reference categories")
        estimator = CoxPHSurvivalAnalysis(
            alpha=config_value(self.config, "alpha"),
            ties=config_value(self.config, "ties"),
            n_iter=config_value(self.config, "n_iter"),
            tol=config_value(self.config, "tol"),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            estimator.fit(active_X, y)
        if not np.isfinite(estimator.coef_).all():
            raise ValueError("Cox fit produced nonfinite coefficients")
        if np.any(np.asarray(config_value(self.config, "alpha")) != 0):
            diagnostic = {
                "status": "unavailable",
                "reason": "ph_score_inference_requires_unpenalized_fit",
                "training_only": True,
            }
        else:
            diagnostic = proportional_hazards_test(
                active_X,
                y,
                estimator.predict(active_X),
                ties=config_value(self.config, "ties"),
                transform=config_value(self.config, "ph_diagnostic_transform"),
                significance_level=config_value(self.config, "ph_diagnostic_alpha"),
            )
        self.estimator_ = estimator
        self.feature_names_ = tuple(X.columns)
        self.active_feature_names_ = active_features
        self.constant_features_ = tuple(name for name in X.columns if name not in active_features)
        self.ph_diagnostic_ = diagnostic
        return self
