"""Predictor transforms learned once on training patients and reused unchanged."""

from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.pipeline import make_pipeline
from sksurv.preprocessing import OneHotEncoder


class TrainingPreprocessor:
    """Median numeric imputation, missingness indicators and named dummy columns.

    Missing and unseen categories use an explicit unknown level. Category sets,
    numeric medians, means and scales are learned exclusively in ``fit``. Numeric
    clinical transforms are identical between the two feature arms.
    """

    def __init__(self, config, feature_set: str):
        self.config = config
        self.feature_set = feature_set
        if feature_set not in {"clinical", "clinical_molecular"}:
            raise ValueError("Unknown feature set.")
        self.features_ = list(config.feature_sets.clinical)
        if feature_set == "clinical_molecular":
            self.features_ += list(config.feature_sets.clinical_molecular)
        self.numeric_ = [
            column for column in self.features_ if column in config.preprocessing.numeric_columns
        ]
        self.categorical_ = [
            column
            for column in self.features_
            if column in config.preprocessing.categorical_columns
        ]
        if set(self.numeric_ + self.categorical_) != set(self.features_):
            raise ValueError("Every predictor must have an explicit preprocessing type.")
        self.unknown_ = config.preprocessing.unknown_category

    @staticmethod
    def _strings(column: pd.Series) -> pd.Series:
        # Grade is ordinal categorical; normalise float CSV values consistently.
        return column.map(
            lambda x: (
                str(int(x))
                if isinstance(x, (float, int)) and pd.notna(x) and float(x).is_integer()
                else (str(x) if pd.notna(x) else None)
            )
        )

    def fit(self, cohort: pd.DataFrame) -> TrainingPreprocessor:
        if self.config.preprocessing.numeric_imputation != "median":
            raise ValueError("The registered preprocessing uses training medians.")
        self.medians_ = {}
        self.means_ = {}
        self.scales_ = {}
        for column in self.numeric_:
            values = cohort[column].astype(float)
            if values.notna().sum() == 0 or np.isinf(values).any():
                raise ValueError("A numeric predictor is entirely missing or contains infinity.")
            self.medians_[column] = float(values.median())
            filled = values.fillna(self.medians_[column])
            self.means_[column] = float(filled.mean())
            scale = float(filled.std(ddof=0))
            self.scales_[column] = scale if scale > 0 else 1.0
        self.categories_ = {}
        for column in self.categorical_:
            observed = sorted(set(self._strings(cohort[column]).dropna()) - {self.unknown_})
            # The first observed level is the dropped reference, so unknown always
            # has a dedicated dummy column when at least one category is observed.
            self.categories_[column] = [*observed, self.unknown_]
        self.encoder_ = make_pipeline(OneHotEncoder())
        encoded = self.encoder_.fit_transform(self._prepare(cohort))
        self.encoded_columns_ = list(encoded.columns)
        return self

    def _prepare(self, cohort: pd.DataFrame) -> pd.DataFrame:
        frame = pd.DataFrame(index=cohort.index)
        for column in self.numeric_:
            values = cohort[column].astype(float)
            frame[column] = (
                values.fillna(self.medians_[column]) - self.means_[column]
            ) / self.scales_[column]
            if self.config.preprocessing.add_missing_indicators:
                frame[f"{column}_missing"] = values.isna().astype(float)
        for column in self.categorical_:
            values = self._strings(cohort[column]).fillna(self.unknown_)
            values = values.where(values.isin(self.categories_[column]), self.unknown_)
            frame[column] = pd.Categorical(values, categories=self.categories_[column])
        return frame

    def transform(self, cohort: pd.DataFrame) -> pd.DataFrame:
        if not hasattr(self, "encoder_"):
            raise ValueError("Preprocessor has not been fitted.")
        output = self.encoder_.transform(self._prepare(cohort)).astype(float)
        if (
            list(output.columns) != self.encoded_columns_
            or not np.isfinite(output.to_numpy()).all()
        ):
            raise ValueError("Transformed feature contract changed or contains non-finite values.")
        return output

    def fit_transform(self, cohort: pd.DataFrame) -> pd.DataFrame:
        return self.fit(cohort).transform(cohort)

    def save(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(pickle.dumps(self, protocol=pickle.HIGHEST_PROTOCOL))

    @classmethod
    def load(cls, path: Path) -> TrainingPreprocessor:
        """Load only this project's locally generated, trusted preprocessing artifact."""
        return pickle.loads(path.read_bytes())
