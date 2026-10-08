"""Shared prediction contract and helpers for classical survival estimators."""

from __future__ import annotations

import pickle
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Protocol, Self

import numpy as np
import pandas as pd
from sklearn.exceptions import NotFittedError


class SurvivalModel(Protocol):
    """All time units are months; higher risk scores mean worse prognosis."""

    def fit(self, X: pd.DataFrame, y: np.ndarray) -> Self: ...
    def predict_risk(self, X: pd.DataFrame) -> np.ndarray: ...
    def predict_survival_function(self, X: pd.DataFrame, times: np.ndarray) -> np.ndarray: ...
    def predict_quantiles(self, X: pd.DataFrame, probabilities: np.ndarray) -> np.ndarray: ...
    def save(self, path: Path | str) -> None: ...

    @classmethod
    def load(cls, path: Path | str) -> Self: ...


def config_value(config: Any, name: str) -> Any:
    """Accept the project's namespace configuration or a plain mapping."""
    return config[name] if isinstance(config, Mapping) else getattr(config, name)


def validate_features(X: pd.DataFrame, expected: tuple[str, ...] | None = None) -> None:
    """Require an already transformed, finite DataFrame with stable named columns."""
    if not isinstance(X, pd.DataFrame) or X.empty:
        raise ValueError("X must be a nonempty pandas DataFrame")
    if not X.columns.is_unique or not all(isinstance(name, str) for name in X.columns):
        raise ValueError("X must have unique string column names")
    if expected is not None and tuple(X.columns) != expected:
        raise ValueError("Prediction columns and order must match the fitted training features")
    if not all(pd.api.types.is_numeric_dtype(dtype) for dtype in X.dtypes):
        raise ValueError("X must contain numeric features after training-only preprocessing")
    if not np.isfinite(X.to_numpy(dtype=float)).all():
        raise ValueError("X must contain only finite, imputed numeric values")


def validate_survival_labels(y: np.ndarray, n_samples: int) -> None:
    """Check the boolean-event / positive floating-point-months Surv contract."""
    if not isinstance(y, np.ndarray) or y.ndim != 1 or len(y) != n_samples:
        raise ValueError("y must be a one-dimensional structured array aligned to X")
    names = y.dtype.names
    if names is None or len(names) != 2:
        raise ValueError("y must have exactly two fields: boolean event and floating-point time")
    if y[names[0]].dtype.kind != "b" or y[names[1]].dtype.kind != "f":
        raise ValueError("y fields must be boolean event followed by floating-point time")
    if not np.isfinite(y[names[1]]).all() or (y[names[1]] <= 0).any():
        raise ValueError("Survival times must be finite and strictly positive months")
    if not y[names[0]].any():
        raise ValueError("Training survival labels must contain at least one observed event")


def step_survival_on_grid(knots: np.ndarray, survival: np.ndarray, times: np.ndarray) -> np.ndarray:
    """Right-continuous steps, 1 before the first knot, last value after the last.

    The carried tail is a computational convention, not identification of survival
    beyond follow-up. Downstream metric grids must stay within training support.
    """
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or not np.isfinite(times).all() or (times < 0).any():
        raise ValueError("Prediction times must be a finite nonnegative one-dimensional array")
    indices = np.searchsorted(knots, times, side="right") - 1
    result = np.ones((survival.shape[0], len(times)), dtype=float)
    known = indices >= 0
    result[:, known] = survival[:, indices[known]]
    return result


def step_survival_quantiles(
    knots: np.ndarray, survival: np.ndarray, probabilities: np.ndarray
) -> np.ndarray:
    """Return inf{t: 1-S(t) >= p}; unreached quantiles are positive infinity.

    Probabilities must lie strictly between zero and one. Infinity records an
    unidentified upper tail; replacing it with the last observed time is invalid.
    """
    probabilities = np.asarray(probabilities, dtype=float)
    if (
        probabilities.ndim != 1
        or not np.isfinite(probabilities).all()
        or ((probabilities <= 0) | (probabilities >= 1)).any()
    ):
        raise ValueError("Quantile probabilities must be a finite vector strictly within (0, 1)")
    result = np.full((survival.shape[0], len(probabilities)), np.inf)
    for j, probability in enumerate(probabilities):
        reached = survival <= 1 - probability
        any_reached = reached.any(axis=1)
        result[any_reached, j] = knots[reached[any_reached].argmax(axis=1)]
    return result


class StepSurvivalModel:
    """Prediction and trusted-artifact persistence for sksurv step estimators."""

    estimator_: Any
    feature_names_: tuple[str, ...]
    active_feature_names_: tuple[str, ...]

    def _check_prediction_input(self, X: pd.DataFrame) -> None:
        if not hasattr(self, "estimator_"):
            raise NotFittedError("Fit the model before predicting or saving")
        validate_features(X, self.feature_names_)

    def predict_risk(self, X: pd.DataFrame) -> np.ndarray:
        self._check_prediction_input(X)
        return np.asarray(
            self.estimator_.predict(X.loc[:, self.active_feature_names_]), dtype=float
        )

    def predict_survival_function(self, X: pd.DataFrame, times: np.ndarray) -> np.ndarray:
        self._check_prediction_input(X)
        survival = self.estimator_.predict_survival_function(
            X.loc[:, self.active_feature_names_], return_array=True
        )
        return step_survival_on_grid(self.estimator_.unique_times_, survival, times)

    def predict_quantiles(self, X: pd.DataFrame, probabilities: np.ndarray) -> np.ndarray:
        self._check_prediction_input(X)
        survival = self.estimator_.predict_survival_function(
            X.loc[:, self.active_feature_names_], return_array=True
        )
        return step_survival_quantiles(self.estimator_.unique_times_, survival, probabilities)

    def save(self, path: Path | str) -> None:
        """Persist a trusted local fit; generated artifacts belong outside version control."""
        if not hasattr(self, "estimator_"):
            raise NotFittedError("Fit the model before saving")
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as stream:
            pickle.dump(
                {"schema_version": 1, "model": self}, stream, protocol=pickle.HIGHEST_PROTOCOL
            )

    @classmethod
    def load(cls, path: Path | str) -> Self:
        """Load only an artifact created by trusted project code (pickle executes code)."""
        with Path(path).open("rb") as stream:
            payload = pickle.load(stream)
        if (
            not isinstance(payload, dict)
            or payload.get("schema_version") != 1
            or type(payload.get("model")) is not cls
        ):
            raise ValueError("Saved model schema or model class does not match")
        return payload["model"]
