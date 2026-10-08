"""Conservative split-conformal lower bounds for a finite restricted event time.

For Z=min(T, horizon), the observed score m(X)-min(Y, horizon) dominates
the latent score m(X)-Z because Y=min(T,C)<=T. The exact finite-sample rank
of the observed calibration scores therefore gives a conservative lower
prediction bound without assuming independent censoring. The model and score
must be fixed before calibration; latent scores must be exchangeable within
each supplied split stratum. Taking the largest stratum threshold avoids
requiring the future patient's unknown event stratum at prediction time.

This is a one-sided score-envelope adaptation, not Qin's bootstrap algorithm
or the published Candes/Gui procedures. The base prediction can be the fitted
model's alpha-quantile; calibration deliberately allows signed corrections.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil
from typing import Any

import numpy as np

from brca.conformal.survival import PredictionIntervals, _labels


def _finite_horizon(horizon: float) -> float:
    if not np.isfinite(horizon) or horizon <= 0:
        raise ValueError("The restricted target requires a finite positive horizon.")
    return float(horizon)


def _restricted_base(base_lower: np.ndarray, horizon: float) -> np.ndarray:
    base = np.asarray(base_lower, dtype=float)
    if base.ndim != 1 or np.any(np.isnan(base)) or np.any(base < 0):
        raise ValueError("Base lower predictions must be a nonnegative vector without NaN.")
    # An unreached model quantile is +infinity for T, but equals the horizon for Z.
    return np.minimum(base, horizon)


def signed_finite_sample_quantile(scores: np.ndarray, alpha: float) -> float:
    """Return the ceil((n+1)(1-alpha))-th signed score, or +infinity."""
    scores = np.asarray(scores, dtype=float)
    if not 0 < alpha < 1:
        raise ValueError("alpha must be miscoverage strictly between zero and one.")
    if scores.ndim != 1 or np.any(~np.isfinite(scores)):
        raise ValueError("Signed calibration scores must be a finite vector.")
    rank = ceil((len(scores) + 1) * (1 - alpha))
    return float(np.sort(scores)[rank - 1]) if rank <= len(scores) else float("inf")


def lower_score_bounds(y_cal: np.ndarray, base_lower: np.ndarray, *, horizon: float) -> np.ndarray:
    """Signed upper bounds m(X)-min(Y,tau) on m(X)-min(T,tau).

    All calibration patients contribute, including censoring before the horizon.
    Values may be negative: truncating scores at zero would needlessly prevent
    calibration from raising an overly conservative base lower prediction.
    """
    horizon = _finite_horizon(horizon)
    if not isinstance(y_cal, np.ndarray) or y_cal.ndim != 1:
        raise ValueError("Calibration outcomes must be a one-dimensional structured array.")
    _, times = _labels(y_cal)
    base = _restricted_base(base_lower, horizon)
    if len(times) != len(base):
        raise ValueError("Outcomes and base lower predictions must be aligned.")
    return base - np.minimum(times, horizon)


@dataclass(frozen=True)
class LowerBoundCalibration:
    alpha: float
    n_cal: int
    horizon: float
    quantile: float
    method: str = "conservative_lower_bound"
    stratum_counts: dict[str, int] = field(default_factory=dict)
    stratum_quantiles: dict[str, float] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _finite_horizon(self.horizon)
        if not 0 < self.alpha < 1 or not isinstance(self.n_cal, int) or self.n_cal < 0:
            raise ValueError("Invalid alpha or calibration-patient count.")
        if np.isnan(self.quantile) or self.quantile == -np.inf:
            raise ValueError("The signed threshold must be finite or positive infinity.")
        if self.method != "conservative_lower_bound":
            raise ValueError("Unexpected method for lower-bound calibration.")

    def predict(self, base_lower: np.ndarray) -> PredictionIntervals:
        """Return [clip(m-q,0,tau), tau]; no patient stratum is required."""
        base = _restricted_base(base_lower, self.horizon)
        if np.isposinf(self.quantile):
            lower = np.zeros_like(base)
        else:
            lower = np.clip(base - self.quantile, 0, self.horizon)
        return PredictionIntervals(
            lower,
            np.full_like(lower, self.horizon),
            self.alpha,
            self.n_cal,
            self.horizon,
            self.method,
            bool(np.isposinf(self.quantile)),
        )

    def to_dict(self) -> dict[str, Any]:
        """JSON-safe calibration state; a null threshold denotes +infinity."""
        return {
            "alpha": self.alpha,
            "n_cal": self.n_cal,
            "horizon": self.horizon,
            "quantile": self.quantile if np.isfinite(self.quantile) else None,
            "quantile_is_infinite": bool(np.isposinf(self.quantile)),
            "method": self.method,
            "stratum_counts": self.stratum_counts,
            "stratum_quantiles": {
                key: value if np.isfinite(value) else None
                for key, value in self.stratum_quantiles.items()
            },
            "diagnostics": self.diagnostics,
        }

    @classmethod
    def from_dict(cls, state: dict[str, Any]) -> LowerBoundCalibration:
        return cls(
            alpha=state["alpha"],
            n_cal=state["n_cal"],
            horizon=state["horizon"],
            quantile=np.inf if state["quantile"] is None else state["quantile"],
            method=state["method"],
            stratum_counts=state.get("stratum_counts", {}),
            stratum_quantiles={
                key: np.inf if value is None else value
                for key, value in state.get("stratum_quantiles", {}).items()
            },
            diagnostics=state.get("diagnostics", {}),
        )


def calibrate_lower_bound(
    y_cal: np.ndarray,
    base_lower: np.ndarray,
    alpha: float,
    *,
    horizon: float,
    strata: np.ndarray | None = None,
    expected_strata: tuple | list | None = None,
) -> LowerBoundCalibration:
    """Calibrate signed lower-bound scores, retaining every calibration patient.

    For the project's event-stratified split, pass primary-endpoint event strata
    for every endpoint. Missing or undersized required strata produce an infinite
    threshold and full restricted support. Neither censoring weights nor future
    outcomes enter the prediction rule.
    """
    horizon = _finite_horizon(horizon)
    scores = lower_score_bounds(y_cal, base_lower, horizon=horizon)
    quantile = signed_finite_sample_quantile(scores, alpha)
    counts, quantiles = {}, {}
    if strata is not None:
        strata = np.asarray(strata).astype(str)
        if strata.shape != scores.shape or expected_strata is None:
            raise ValueError("Aligned strata and explicit expected_strata are required together.")
        levels = [str(level) for level in expected_strata]
        if not levels or len(set(levels)) != len(levels) or set(strata) - set(levels):
            raise ValueError("expected_strata must contain every allowed stratum exactly once.")
        for level in levels:
            selected = scores[strata == level]
            counts[level] = len(selected)
            quantiles[level] = signed_finite_sample_quantile(selected, alpha)
        quantile = max(quantiles.values())
    elif expected_strata is not None:
        raise ValueError("expected_strata requires calibration strata.")
    return LowerBoundCalibration(
        alpha=alpha,
        n_cal=len(scores),
        horizon=horizon,
        quantile=quantile,
        stratum_counts=counts,
        stratum_quantiles=quantiles,
        diagnostics={
            "n_censored_before_horizon": int((~y_cal["event"] & (y_cal["time"] < horizon)).sum()),
            "n_negative_scores": int((scores < 0).sum()),
            "score": "min(base_lower,horizon)-min(observed_time,horizon)",
            "target": "min(event_time,horizon)",
            "interval_type": "one_sided_lower_bound_with_fixed_horizon_upper_edge",
            "censoring_assumption": "none_for_pathwise_score_domination",
            "rank_assumption": "exchangeability_within_split_strata"
            if strata is not None
            else "exchangeability",
        },
    )
