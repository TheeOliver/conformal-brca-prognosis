"""Untuned random survival forest, with an injected reproducible random stream."""

from __future__ import annotations

from typing import Any, Self

import numpy as np
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest

from brca.models.base import (
    StepSurvivalModel,
    config_value,
    validate_features,
    validate_survival_labels,
)


class RandomSurvivalForestModel(StepSurvivalModel):
    """Fixed configured hyperparameters; no calibration/test selection or tuning."""

    def __init__(self, config: Any, *, rng: np.random.Generator):
        self.config = config
        self.random_state = int(rng.integers(0, np.iinfo(np.int32).max))

    def fit(self, X: pd.DataFrame, y: np.ndarray) -> Self:
        validate_features(X)
        validate_survival_labels(y, len(X))
        n_jobs = config_value(self.config, "n_jobs")
        if not isinstance(n_jobs, int) or n_jobs < 1:
            raise ValueError("RSF n_jobs must be a positive allocated CPU count, never -1")
        estimator = RandomSurvivalForest(
            n_estimators=config_value(self.config, "n_estimators"),
            min_samples_split=config_value(self.config, "min_samples_split"),
            min_samples_leaf=config_value(self.config, "min_samples_leaf"),
            max_features=config_value(self.config, "max_features"),
            n_jobs=n_jobs,
            random_state=self.random_state,
        )
        estimator.fit(X, y)
        self.estimator_ = estimator
        self.feature_names_ = tuple(X.columns)
        self.active_feature_names_ = tuple(X.columns)
        return self
