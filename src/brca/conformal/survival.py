"""Conservative split conformal intervals for right-censored survival outcomes.

The primary target is Z=min(T, horizon). Calibration uses the supremum of a
nonnegative CQR score over every survival time compatible with each observation.
Its rank argument requires exchangeability within the specified split strata,
but does not require independent censoring. Stratum thresholds are maximised so
prediction never requires knowing the future patient's event indicator.

IPCW calibration and coverage estimates are separate sensitivity analyses. Their
marginal reverse-KM weights assume C independent of (T, X), positive support and
consistent censoring estimation; they have no asserted finite-sample guarantee.
See docs/conformal-methods.md for the estimand, proof and competing-risk caveat.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil
from typing import Any

import numpy as np
from scipy.stats import norm
from sksurv.nonparametric import kaplan_meier_estimator


def _labels(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if y.dtype.names != ("event", "time") or y["event"].dtype.kind != "b":
        raise ValueError("y must have structured fields ('event', 'time'), with bool event.")
    time = np.asarray(y["time"], dtype=float)
    if time.ndim != 1 or np.any(~np.isfinite(time)) or np.any(time <= 0):
        raise ValueError("Survival times must be finite, strictly positive months.")
    return y["event"], time


def _horizon(horizon: float) -> float:
    if np.isnan(horizon) or horizon <= 0:
        raise ValueError("horizon must be positive; infinity denotes unrestricted survival.")
    return float(horizon)


def _bounds(lower: np.ndarray, upper: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    if lower.ndim != 1 or lower.shape != upper.shape:
        raise ValueError("lower and upper must be aligned one-dimensional arrays.")
    if np.any(~np.isfinite(lower)) or np.any(np.isnan(upper)):
        raise ValueError("Lower endpoints must be finite and upper endpoints cannot be NaN.")
    if np.any(lower < 0) or np.any(upper < lower):
        raise ValueError("Intervals must satisfy 0 <= lower <= upper.")
    return lower, upper


def finite_sample_quantile(scores: np.ndarray, alpha: float) -> float:
    """Return the ceil((n+1)(1-alpha))-th score; missing rank means +infinity."""
    scores = np.asarray(scores, float)
    if not 0 < alpha < 1:
        raise ValueError("alpha is miscoverage and must be between zero and one.")
    if scores.ndim != 1 or np.any(np.isnan(scores)) or np.any(scores < 0):
        raise ValueError("scores must be a vector of nonnegative, non-NaN values.")
    k = ceil((len(scores) + 1) * (1 - alpha))
    return float(np.sort(scores)[k - 1]) if k <= len(scores) else float("inf")


def score_bounds(
    y: np.ndarray, lower: np.ndarray, upper: np.ndarray, *, horizon: float
) -> np.ndarray:
    """Upper-bound the latent nonnegative CQR score, retaining censored patients.

    Known Z has max(lower-Z, Z-upper, 0). If censoring precedes the horizon,
    the largest compatible score is max(lower-observed_time, horizon-upper, 0).
    Unrestricted censored observations conservatively receive infinity.
    """
    event, time = _labels(y)
    horizon = _horizon(horizon)
    lower, upper = _bounds(lower, upper)
    if len(time) != len(lower):
        raise ValueError("Outcomes and prediction bounds must have the same length.")
    lower, upper = np.minimum(lower, horizon), np.minimum(upper, horizon)
    z = np.minimum(time, horizon)
    known = event | (time >= horizon)
    scores = np.maximum(np.maximum(lower - z, z - upper), 0)
    if np.isfinite(horizon):
        scores[~known] = np.maximum.reduce(
            [lower[~known] - z[~known], horizon - upper[~known], np.zeros((~known).sum())]
        )
    else:
        scores[~known] = np.inf
    return scores


@dataclass(frozen=True)
class PredictionIntervals:
    lower: np.ndarray
    upper: np.ndarray
    alpha: float
    n_cal: int
    horizon: float
    method: str = "model_predictive"
    quantile_is_infinite: bool = False

    def __post_init__(self) -> None:
        lower, upper = _bounds(self.lower, self.upper)
        _horizon(self.horizon)
        if not 0 < self.alpha < 1 or self.n_cal < 0:
            raise ValueError("Invalid alpha or calibration count.")
        if np.any(upper > self.horizon):
            raise ValueError("Interval endpoints must respect the target horizon.")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @property
    def width(self) -> np.ndarray:
        return self.upper - self.lower

    def subset(self, mask: np.ndarray) -> PredictionIntervals:
        return PredictionIntervals(
            self.lower[mask],
            self.upper[mask],
            self.alpha,
            self.n_cal,
            self.horizon,
            self.method,
            self.quantile_is_infinite,
        )


@dataclass(frozen=True)
class ConformalCalibration:
    alpha: float
    n_cal: int
    horizon: float
    quantile: float
    method: str
    stratum_counts: dict[str, int] = field(default_factory=dict)
    stratum_quantiles: dict[str, float] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def predict(self, lower: np.ndarray, upper: np.ndarray) -> PredictionIntervals:
        lower, upper = _bounds(lower, upper)
        if np.isinf(self.quantile):
            lo, hi = np.zeros_like(lower), np.full_like(upper, self.horizon)
        else:
            lo = np.maximum(0, np.minimum(lower, self.horizon) - self.quantile)
            hi = np.minimum(self.horizon, np.minimum(upper, self.horizon) + self.quantile)
        return PredictionIntervals(
            lo,
            hi,
            self.alpha,
            self.n_cal,
            self.horizon,
            self.method,
            bool(np.isinf(self.quantile)),
        )

    def to_dict(self) -> dict[str, Any]:
        """JSON-safe state; null thresholds/horizon mean positive infinity."""
        return {
            "alpha": self.alpha,
            "n_cal": self.n_cal,
            "horizon": self.horizon if np.isfinite(self.horizon) else None,
            "quantile": self.quantile if np.isfinite(self.quantile) else None,
            "quantile_is_infinite": bool(np.isinf(self.quantile)),
            "method": self.method,
            "stratum_counts": self.stratum_counts,
            "stratum_quantiles": {
                key: value if np.isfinite(value) else None
                for key, value in self.stratum_quantiles.items()
            },
            "diagnostics": self.diagnostics,
        }

    @classmethod
    def from_dict(cls, state: dict[str, Any]) -> ConformalCalibration:
        return cls(
            alpha=state["alpha"],
            n_cal=state["n_cal"],
            horizon=np.inf if state["horizon"] is None else state["horizon"],
            quantile=np.inf if state["quantile"] is None else state["quantile"],
            method=state["method"],
            stratum_counts=state.get("stratum_counts", {}),
            stratum_quantiles={
                key: np.inf if value is None else value
                for key, value in state.get("stratum_quantiles", {}).items()
            },
            diagnostics=state.get("diagnostics", {}),
        )


def calibrate_intervals(
    y_cal: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    alpha: float,
    *,
    horizon: float,
    strata: np.ndarray | None = None,
    expected_strata: tuple | list | None = None,
) -> ConformalCalibration:
    """Calibrate conservative intervals; fix horizon, score and strata before calibration.

    Supply the primary-event stratum even for a secondary endpoint. Maximise the
    within-stratum quantiles to avoid conditioning prediction on an unknown label.
    A missing required stratum contributes infinity, never a fabricated threshold.
    """
    scores = score_bounds(y_cal, lower, upper, horizon=horizon)
    quantile = finite_sample_quantile(scores, alpha)
    counts, quantiles = {}, {}
    if strata is not None:
        strata = np.asarray(strata).astype(str)
        if strata.shape != scores.shape or expected_strata is None:
            raise ValueError("Aligned strata and explicit expected_strata are required together.")
        levels = [str(level) for level in expected_strata]
        if not levels or len(set(levels)) != len(levels) or set(strata) - set(levels):
            raise ValueError("expected_strata must list every allowed stratum exactly once.")
        for level in levels:
            selected = scores[strata == level]
            counts[level] = len(selected)
            quantiles[level] = finite_sample_quantile(selected, alpha)
        quantile = max(quantiles.values())
    elif expected_strata is not None:
        raise ValueError("expected_strata requires calibration strata.")
    return ConformalCalibration(
        alpha,
        len(scores),
        horizon,
        quantile,
        "conservative_score_envelope",
        counts,
        quantiles,
        {
            "n_censored_before_horizon": int((~y_cal["event"] & (y_cal["time"] < horizon)).sum()),
            "rank_assumption": "exchangeability_within_split_strata"
            if strata is not None
            else "exchangeability",
        },
    )


@dataclass(frozen=True)
class CensoringDistribution:
    """Reverse Kaplan-Meier fit on train only, with event-before-censor tie order."""

    times: np.ndarray
    survival_after: np.ndarray

    @classmethod
    def fit(cls, y_train: np.ndarray) -> CensoringDistribution:
        event, time = _labels(y_train)
        if not len(time):
            raise ValueError("Training outcomes are needed to estimate censoring.")
        times, survival = kaplan_meier_estimator(event, time, reverse=True)
        return cls(times, survival)

    def survival_ge(self, time: np.ndarray) -> np.ndarray:
        """Estimate P(C >= t)=G(t-) without silently extrapolating past train support."""
        time = np.asarray(time, float)
        index = np.searchsorted(self.times, time, side="left") - 1
        values = np.ones_like(time, dtype=float)
        present = index >= 0
        values[present] = self.survival_after[index[present]]
        values[time > self.times[-1]] = np.nan
        return values


def _ipcw_weights(
    y: np.ndarray,
    censoring: CensoringDistribution,
    horizon: float,
    min_censoring_survival: float,
) -> tuple[np.ndarray | None, dict[str, Any]]:
    event, time = _labels(y)
    if not 0 < min_censoring_survival <= 1:
        raise ValueError("min_censoring_survival must be in (0, 1].")
    known = event | (time >= horizon)
    diagnostics = {"n_known_target": int(known.sum()), "status": "available"}
    if not np.isfinite(horizon):
        diagnostics["status"] = "unavailable_unrestricted_tail"
        return None, diagnostics
    g_horizon = float(censoring.survival_ge(np.array([horizon]))[0])
    diagnostics["censoring_survival_at_horizon"] = g_horizon if np.isfinite(g_horizon) else None
    if not np.isfinite(g_horizon) or g_horizon < min_censoring_survival:
        diagnostics["status"] = "unavailable_censoring_support"
        return None, diagnostics
    weights = np.zeros(len(time))
    weights[known] = 1 / censoring.survival_ge(np.minimum(time[known], horizon))
    diagnostics["max_weight"] = float(weights.max(initial=0))
    diagnostics["effective_sample_size"] = (
        float(weights.sum() ** 2 / np.square(weights).sum()) if weights.sum() else 0.0
    )
    return weights, diagnostics


def calibrate_ipcw_intervals(
    y_cal: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    y_train: np.ndarray,
    alpha: float,
    *,
    horizon: float,
    min_censoring_survival: float,
) -> ConformalCalibration:
    """Optional approximate IPCW sensitivity, not the primary conformal guarantee.

    Use observed restricted targets with inverse censoring weights, plus an
    infinity atom of weight 1/G(horizon-). In the fully observed case every
    weight is one, reducing exactly to the split conformal order statistic.
    Estimated weights, event-stratified sampling and marginal rather than
    conditional censoring estimation preclude asserting finite-sample validity.
    """
    scores = score_bounds(y_cal, lower, upper, horizon=horizon)
    finite_sample_quantile(scores, alpha)  # validate alpha even if support is absent
    weights, diagnostics = _ipcw_weights(
        y_cal, CensoringDistribution.fit(y_train), horizon, min_censoring_survival
    )
    diagnostics["validity"] = "approximate_sensitivity_no_finite_sample_guarantee"
    quantile = np.inf
    if weights is not None and weights.sum() > 0:
        included = weights > 0
        order = np.argsort(scores[included], kind="stable")
        ordered_scores, ordered_weights = scores[included][order], weights[included][order]
        infinity_weight = 1 / diagnostics["censoring_survival_at_horizon"]
        target = (1 - alpha) * (weights.sum() + infinity_weight)
        position = np.searchsorted(np.cumsum(ordered_weights), target, side="left")
        if position < len(ordered_scores):
            quantile = float(ordered_scores[position])
        diagnostics["infinity_atom_weight"] = infinity_weight
    return ConformalCalibration(
        alpha, len(scores), horizon, quantile, "ipcw_score_sensitivity", diagnostics=diagnostics
    )


def _wilson(successes: int, n: int, confidence_level: float) -> list[float] | None:
    if not n:
        return None
    z = norm.ppf((1 + confidence_level) / 2)
    p, correction = successes / n, z * z / n
    center = (p + correction / 2) / (1 + correction)
    radius = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + correction)
    return [float(max(0, center - radius)), float(min(1, center + radius))]


def _finite_number(value: float) -> float | None:
    return float(value) if np.isfinite(value) else None


def _coverage_indicators(
    y: np.ndarray, intervals: PredictionIntervals
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    event, time = _labels(y)
    if len(time) != len(intervals.lower):
        raise ValueError("Test labels and prediction intervals must be aligned.")
    z = np.minimum(time, intervals.horizon)
    known = event | (time >= intervals.horizon)
    covered = (intervals.lower <= z) & (z <= intervals.upper)
    # Censoring implies T > Y; its compatible restricted values are (Y, tau].
    certain = np.where(
        known, covered, (intervals.lower <= z) & (intervals.upper >= intervals.horizon)
    )
    possible = np.where(
        known, covered, (intervals.upper > z) & (intervals.lower <= intervals.horizon)
    )
    return certain, possible, covered


def _bootstrap_interval_draws(
    y: np.ndarray,
    intervals: PredictionIntervals,
    selections: np.ndarray,
    weights: np.ndarray | None,
) -> dict[str, np.ndarray]:
    """Use precisely the same patient selections for every reported statistic."""
    certain, possible, covered = _coverage_indicators(y, intervals)
    width_quantiles = np.quantile(
        intervals.width[selections], [0.25, 0.5, 0.75], axis=1, method="inverted_cdf"
    )
    full_support = (intervals.lower == 0) & (intervals.upper == intervals.horizon)
    draws = {
        "observable_coverage_lower": certain[selections].mean(axis=1),
        "observable_coverage_upper": possible[selections].mean(axis=1),
        "width_q25_months": width_quantiles[0],
        "width_median_months": width_quantiles[1],
        "width_q75_months": width_quantiles[2],
        "full_support_fraction": full_support[selections].mean(axis=1),
        "infinite_interval_fraction": np.isinf(intervals.width)[selections].mean(axis=1),
        "ipcw_coverage": np.full(len(selections), np.nan),
        "ipcw_ht_coverage": np.full(len(selections), np.nan),
    }
    if weights is not None:
        denominator = weights[selections].sum(axis=1)
        numerator = (weights * covered)[selections].sum(axis=1)
        np.divide(numerator, denominator, out=draws["ipcw_coverage"], where=denominator > 0)
        draws["ipcw_ht_coverage"] = numerator / len(y)
    return draws


def _validate_bootstrap(
    replicates: int,
    confidence_level: float,
    minimum_success_fraction: float,
    strata: np.ndarray | None,
    n: int,
) -> np.ndarray | None:
    if (
        not isinstance(replicates, int)
        or isinstance(replicates, bool)
        or replicates < 0
        or not 0 < confidence_level < 1
        or not 0 < minimum_success_fraction <= 1
    ):
        raise ValueError("Invalid bootstrap count, confidence level or minimum success fraction.")
    if strata is not None:
        strata = np.asarray(strata)
        if strata.shape != (n,):
            raise ValueError("bootstrap_strata must align with test observations.")
    return strata


def _bootstrap_indices(
    n: int,
    replicates: int,
    seed: int,
    strata: np.ndarray | None,
) -> np.ndarray:
    """Resample patients within their original split strata, keeping stratum counts fixed."""
    rng = np.random.default_rng(seed)
    if strata is None:
        return rng.integers(0, n, size=(replicates, n))
    labels = np.asarray(strata).astype(str)
    selections = np.empty((replicates, n), dtype=int)
    for level in np.unique(labels):
        members = np.flatnonzero(labels == level)
        selections[:, members] = rng.choice(members, size=(replicates, len(members)))
    return selections


def _draw_summary(
    draws: np.ndarray,
    confidence_level: float,
    minimum_success_fraction: float,
) -> dict[str, Any]:
    # Keep infinities explicit; they are not missing estimates. Inverted-CDF
    # quantiles avoid interpolation of infinite widths into NaNs.
    requested = len(draws)
    draws = draws[~np.isnan(draws)]
    n_valid = len(draws)
    fraction = n_valid / requested if requested else None
    passed = n_valid >= 2 and fraction >= minimum_success_fraction
    result = {
        "ci": None,
        "ci_is_infinite": None,
        "valid_replicates": n_valid,
        "failed_replicates": requested - n_valid,
        "success_fraction": fraction,
        "minimum_success_fraction": minimum_success_fraction,
        "status": "bootstrap_gate_failed" if requested else "not_requested",
        "approximate_due_to_failed_resamples": False,
    }
    if not passed:
        return result
    tail = (1 - confidence_level) / 2
    limits = np.quantile(draws, [tail, 1 - tail], method="inverted_cdf")
    result.update(
        {
            "ci": [_finite_number(value) for value in limits],
            "ci_is_infinite": np.isinf(limits).tolist(),
            "status": "ok" if n_valid == requested else "conditional_on_successful_resamples",
            "approximate_due_to_failed_resamples": n_valid < requested,
        }
    )
    return result


def coverage_summary(
    y_test: np.ndarray,
    intervals: PredictionIntervals,
    y_train: np.ndarray,
    *,
    bootstrap_replicates: int,
    seed: int,
    confidence_level: float,
    min_censoring_survival: float,
    minimum_bootstrap_success_fraction: float,
    bootstrap_strata: np.ndarray | None = None,
) -> dict[str, Any]:
    """Summarise coverage bounds, assumption-dependent IPCW coverage and width.

    Bounds contain the realised latent test coverage without independent
    censoring. Wilson intervals describe the two observable Bernoulli rates,
    not a binomial CI for latent coverage. IPCW percentile intervals bootstrap
    test patients conditional on the fitted train censoring distribution, base
    model and calibration state. They omit nuisance/model/refitting uncertainty.
    The same seed and aligned rows permit paired test resampling across arms.
    Supply primary DSS event strata for the fixed event-stratified test split,
    including for OS. CIs with discarded draws are explicitly approximate and
    conditional on successful resamples, subject to the success-fraction gate.
    """
    event, time = _labels(y_test)
    n = len(time)
    if n != len(intervals.lower):
        raise ValueError("Test labels and prediction intervals must be aligned.")
    bootstrap_strata = _validate_bootstrap(
        bootstrap_replicates,
        confidence_level,
        minimum_bootstrap_success_fraction,
        bootstrap_strata,
        n,
    )
    horizon = intervals.horizon
    known = event | (time >= horizon)
    certain, possible, covered = _coverage_indicators(y_test, intervals)
    widths = intervals.width
    # Order statistics avoid invalid infinity interpolation (inf-inf -> NaN).
    width_quantiles = np.quantile(widths, [0.25, 0.5, 0.75], method="inverted_cdf") if n else None
    report: dict[str, Any] = {
        "method": intervals.method,
        "alpha": intervals.alpha,
        "nominal_coverage": 1 - intervals.alpha,
        "n_cal": intervals.n_cal,
        "n_test": n,
        "horizon_months": _finite_number(horizon),
        "target": "min(T,horizon)" if np.isfinite(horizon) else "T",
        "n_known_target": int(known.sum()),
        "n_censored_before_horizon": int((~known).sum()),
        "observable_coverage_lower": float(certain.mean()) if n else None,
        "observable_coverage_upper": float(possible.mean()) if n else None,
        "observable_lower_wilson_ci": _wilson(int(certain.sum()), n, confidence_level),
        "observable_upper_wilson_ci": _wilson(int(possible.sum()), n, confidence_level),
        "width_median_months": _finite_number(width_quantiles[1]) if n else None,
        "width_q25_months": _finite_number(width_quantiles[0]) if n else None,
        "width_q75_months": _finite_number(width_quantiles[2]) if n else None,
        "width_median_is_infinite": bool(n and np.isinf(width_quantiles[1])),
        "width_q25_is_infinite": bool(n and np.isinf(width_quantiles[0])),
        "width_q75_is_infinite": bool(n and np.isinf(width_quantiles[2])),
        "infinite_interval_fraction": float(np.isinf(widths).mean()) if n else None,
        "full_support_fraction": float(
            ((intervals.lower == 0) & (intervals.upper == horizon)).mean()
        )
        if n
        else None,
        "quantile_is_infinite": intervals.quantile_is_infinite,
        "ipcw_coverage": None,
        "ipcw_coverage_ci": None,
        "ipcw_ht_coverage": None,
        "ipcw_ht_coverage_ci": None,
        "observable_coverage_lower_ci": None,
        "observable_coverage_upper_ci": None,
        "width_median_ci_months": None,
        "width_q25_ci_months": None,
        "width_q75_ci_months": None,
        "full_support_fraction_ci": None,
        "infinite_interval_fraction_ci": None,
        "bootstrap_statistics": {},
        "confidence_level": confidence_level,
        "bootstrap_replicates_requested": bootstrap_replicates,
        "bootstrap_replicates_valid": 0,
        "minimum_bootstrap_success_fraction": minimum_bootstrap_success_fraction,
        "bootstrap_sampling": "within_primary_event_strata"
        if bootstrap_strata is not None
        else "ordinary_test_patient",
        "ci_scope": "test_bootstrap_conditional_on_model_calibration_and_train_censoring_fit",
        "ipcw_assumption": "C independent of (T,X), positive censoring support",
    }
    if not n:
        report["ipcw_status"] = "unavailable_empty_subgroup"
        return report
    weights, diagnostics = _ipcw_weights(
        y_test, CensoringDistribution.fit(y_train), horizon, min_censoring_survival
    )
    report["ipcw_status"] = diagnostics.pop("status")
    report["ipcw_diagnostics"] = diagnostics
    if weights is not None and weights.sum() == 0:
        report["ipcw_status"] = "unavailable_no_observed_targets"
        weights = None
    if weights is not None:
        contributions = weights * covered
        report["ipcw_coverage"] = float(contributions.sum() / weights.sum())
        report["ipcw_ht_coverage"] = float(contributions.mean())
        report["ipcw_observed_mass"] = float(weights.mean())
    if bootstrap_replicates:
        selections = _bootstrap_indices(n, bootstrap_replicates, seed, bootstrap_strata)
        draws = _bootstrap_interval_draws(y_test, intervals, selections, weights)
        summaries = {
            key: _draw_summary(value, confidence_level, minimum_bootstrap_success_fraction)
            for key, value in draws.items()
        }
        if weights is None:
            for key in ("ipcw_coverage", "ipcw_ht_coverage"):
                summaries[key]["status"] = report["ipcw_status"]
        report["bootstrap_statistics"] = summaries
        report["bootstrap_replicates_valid"] = summaries["ipcw_coverage"]["valid_replicates"]
        for key, summary in summaries.items():
            ci_key = (
                key.replace("_months", "_ci_months") if key.endswith("_months") else f"{key}_ci"
            )
            report[ci_key] = summary["ci"]
    return report


def paired_interval_comparison(
    y_test: np.ndarray,
    intervals_clinical: PredictionIntervals,
    intervals_molecular: PredictionIntervals,
    y_train: np.ndarray,
    *,
    bootstrap_replicates: int,
    seed: int,
    confidence_level: float,
    min_censoring_survival: float,
    minimum_bootstrap_success_fraction: float,
    bootstrap_strata: np.ndarray | None = None,
) -> dict[str, Any]:
    """Molecular-minus-clinical differences with paired test-patient bootstrap.

    Both arms must refer to the same ordered patients, alpha and target horizon.
    Width compares the arms' marginal medians (not the median of individual
    width differences); positive means wider molecular intervals. Censoring
    weights are shared and fit on train. CIs condition on all fitted objects.
    """
    event, _time = _labels(y_test)
    n = len(event)
    if (
        intervals_clinical.alpha != intervals_molecular.alpha
        or intervals_clinical.horizon != intervals_molecular.horizon
    ):
        raise ValueError("Paired comparisons require the same alpha and horizon.")
    if n != len(intervals_clinical.lower) or n != len(intervals_molecular.lower):
        raise ValueError("Both interval arms and test outcomes must be aligned.")
    bootstrap_strata = _validate_bootstrap(
        bootstrap_replicates,
        confidence_level,
        minimum_bootstrap_success_fraction,
        bootstrap_strata,
        n,
    )
    horizon = intervals_clinical.horizon
    report: dict[str, Any] = {
        "direction": "clinical_molecular_minus_clinical",
        "n_test": n,
        "n_cal_clinical": intervals_clinical.n_cal,
        "n_cal_molecular": intervals_molecular.n_cal,
        "alpha": intervals_clinical.alpha,
        "horizon_months": _finite_number(horizon),
        "same_patient_resamples": True,
        "confidence_level": confidence_level,
        "bootstrap_replicates_requested": bootstrap_replicates,
        "minimum_bootstrap_success_fraction": minimum_bootstrap_success_fraction,
        "bootstrap_sampling": "within_primary_event_strata"
        if bootstrap_strata is not None
        else "ordinary_test_patient",
        "ci_scope": "test_bootstrap_conditional_on_model_calibration_and_train_censoring_fit",
        "width_interpretation": "positive_difference_means_wider_molecular_intervals",
        "differences": {},
    }
    keys = [
        "ipcw_coverage",
        "observable_coverage_lower",
        "observable_coverage_upper",
        "width_median_months",
        "width_q25_months",
        "width_q75_months",
        "full_support_fraction",
        "infinite_interval_fraction",
        "ipcw_ht_coverage",
    ]
    if not n:
        report["ipcw_status"] = "unavailable_empty_subgroup"
        report["differences"] = {
            key: {
                "estimate": None,
                "ci": None,
                "valid_replicates": 0,
                "failed_replicates": 0,
                "success_fraction": None,
                "status": "unavailable_empty_subgroup",
            }
            for key in keys
        }
        return report
    weights, diagnostics = _ipcw_weights(
        y_test, CensoringDistribution.fit(y_train), horizon, min_censoring_survival
    )
    report["ipcw_status"] = diagnostics.pop("status")
    if weights is not None and weights.sum() == 0:
        weights = None
        report["ipcw_status"] = "unavailable_no_observed_targets"
    report["ipcw_diagnostics"] = diagnostics
    # Row zero evaluates the point statistic without resampling. Remaining rows
    # use identical selections for the two arms and every statistic.
    selections = np.vstack(
        [np.arange(n), _bootstrap_indices(n, bootstrap_replicates, seed, bootstrap_strata)]
    )
    clinical = _bootstrap_interval_draws(y_test, intervals_clinical, selections, weights)
    molecular = _bootstrap_interval_draws(y_test, intervals_molecular, selections, weights)
    for key in keys:
        # infinity-infinity is genuinely undefined, not a zero width difference.
        with np.errstate(invalid="ignore"):
            differences = molecular[key] - clinical[key]
        summary = _draw_summary(
            differences[1:], confidence_level, minimum_bootstrap_success_fraction
        )
        if weights is None and key in ("ipcw_coverage", "ipcw_ht_coverage"):
            summary["status"] = report["ipcw_status"]
        summary["estimate"] = _finite_number(differences[0])
        summary["estimate_is_infinite"] = bool(np.isinf(differences[0]))
        summary["estimate_is_undefined"] = bool(np.isnan(differences[0]))
        if np.isinf(differences[0]):
            summary["infinite_direction"] = "positive" if differences[0] > 0 else "negative"
        report["differences"][key] = summary
    return report


def subgroup_coverage(
    y_test: np.ndarray,
    intervals: PredictionIntervals,
    y_train: np.ndarray,
    *,
    groups: dict[str, np.ndarray],
    expected_levels: dict[str, list] | None = None,
    bootstrap_strata: np.ndarray | None = None,
    **summary_kwargs: Any,
) -> dict[str, dict[str, dict[str, Any]]]:
    """Use the same marginal intervals in every subgroup; do not recalibrate.

    Explicit missing levels receive n_test=0 and null estimates. Labels missing
    in a subgroup column become the literal 'unknown'. Small subgroup estimates
    remain visible with their sample count and uncertainty, without validity claims.
    """
    if bootstrap_strata is not None:
        bootstrap_strata = np.asarray(bootstrap_strata)
        if bootstrap_strata.shape != (len(y_test),):
            raise ValueError("bootstrap_strata must align with test observations.")
    result = {}
    for name, values in groups.items():
        values = np.asarray(values)
        if values.shape != (len(y_test),):
            raise ValueError("Subgroup labels must align with the test observations.")
        labels = np.array(
            [
                "unknown"
                if value is None or (isinstance(value, float) and np.isnan(value))
                else str(value)
                for value in values
            ]
        )
        levels = list(
            dict.fromkeys(
                [*(str(x) for x in (expected_levels or {}).get(name, [])), *sorted(set(labels))]
            )
        )
        result[name] = {}
        for level in levels:
            mask = labels == level
            result[name][level] = coverage_summary(
                y_test[mask],
                intervals.subset(mask),
                y_train,
                bootstrap_strata=bootstrap_strata[mask] if bootstrap_strata is not None else None,
                **summary_kwargs,
            )
    return result
