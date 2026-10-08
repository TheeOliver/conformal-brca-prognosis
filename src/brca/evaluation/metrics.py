"""Training-supported IPCW metrics and paired test-patient bootstrap intervals."""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from contextlib import suppress
from typing import Any

import numpy as np
from sksurv.exceptions import NoComparablePairException
from sksurv.metrics import (
    brier_score,
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    integrated_brier_score,
)
from sksurv.nonparametric import CensoringDistributionEstimator

from brca.models.base import config_value, validate_survival_labels

SCALAR_METRICS = ("uno_c", "ibs", "mean_auc")
CURVE_METRICS = ("brier", "auc")


class EvaluationSupportError(ValueError):
    """The prespecified evaluation cannot be estimated on the held-out support."""

    def __init__(self, message: str, diagnostics: dict[str, Any]):
        super().__init__(message)
        self.diagnostics = diagnostics


def _settings(config: Any) -> Any:
    if isinstance(config, Mapping):
        return config.get("evaluation", config)
    return getattr(config, "evaluation", config)


def training_time_grid(y_train: np.ndarray, config: Any) -> np.ndarray:
    """Build one grid from training follow-up alone, capped by the protocol horizon.

    Percentiles above the horizon are capped, then duplicate knots are removed.
    The last knot must precede the final training event; no held-out information
    is available here to choose or move the grid.
    """
    validate_survival_labels(y_train, len(y_train))
    settings = _settings(config)
    percentiles = np.asarray(config_value(settings, "time_grid_percentiles"), dtype=float)
    horizon = float(config_value(settings, "maximum_horizon_months"))
    if (
        percentiles.ndim != 1
        or len(percentiles) < 2
        or not np.isfinite(percentiles).all()
        or ((percentiles <= 0) | (percentiles >= 100)).any()
        or not np.isfinite(horizon)
        or horizon <= 0
    ):
        raise ValueError(
            "Evaluation needs at least two interior percentiles and a positive horizon"
        )
    event, time = (y_train[name] for name in y_train.dtype.names)
    upper = min(horizon, np.nextafter(time[event].max(), -np.inf))
    times = np.unique(np.minimum(np.percentile(time, percentiles), upper))
    times = times[times > time.min()]
    if len(times) < 2:
        raise EvaluationSupportError(
            "Training data cannot support two distinct interior evaluation times",
            {"n_training": len(y_train), "n_training_events": int(event.sum())},
        )
    censoring = CensoringDistributionEstimator().fit(y_train)
    if (censoring.predict_proba(np.append(times, np.nextafter(times[-1], np.inf))) <= 0).any():
        raise EvaluationSupportError(
            "Training censoring survival is zero on the proposed evaluation support",
            {"time_grid_months": times.tolist()},
        )
    return times


def _capped_test_labels(y_test: np.ndarray, horizon: float) -> tuple[np.ndarray, int]:
    """Preserve all cases through horizon and retain later patients as controls."""
    capped = y_test.copy()
    event_name, time_name = y_test.dtype.names
    cap = np.nextafter(horizon, np.inf)
    after_horizon = capped[time_name] > horizon
    capped[event_name][after_horizon] = False
    capped[time_name][after_horizon] = cap
    return capped, int(after_horizon.sum())


def _metric_values(
    y_train: np.ndarray,
    y_test: np.ndarray,
    risk: np.ndarray,
    survival: np.ndarray,
    times: np.ndarray,
) -> dict[str, np.ndarray]:
    """Undefined bootstrap metrics remain NaN internally and are counted downstream."""
    values = {name: np.array(np.nan) for name in SCALAR_METRICS}
    values.update({name: np.full(len(times), np.nan) for name in CURVE_METRICS})
    exceptions = (ValueError, NoComparablePairException, ZeroDivisionError, RuntimeWarning)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with suppress(*exceptions):
            values["uno_c"] = np.array(
                concordance_index_ipcw(y_train, y_test, risk, tau=np.nextafter(times[-1], np.inf))[
                    0
                ]
            )
        with suppress(*exceptions):
            values["brier"] = brier_score(y_train, y_test, survival, times)[1]
            values["ibs"] = np.array(integrated_brier_score(y_train, y_test, survival, times))
        with suppress(*exceptions):
            # 1-S(t) is an event risk at t. Supplying S(t) would reverse AUC ordering.
            values["auc"], mean_auc = cumulative_dynamic_auc(y_train, y_test, 1 - survival, times)
            values["mean_auc"] = np.array(mean_auc)
    return values


def _interval_summary(
    estimate: float,
    replicates: np.ndarray,
    confidence_level: float,
    minimum_success_fraction: float,
) -> dict[str, Any]:
    valid = np.isfinite(replicates)
    n_valid = int(valid.sum())
    fraction = n_valid / len(replicates)
    passed = fraction >= minimum_success_fraction and n_valid >= 2
    bounds = [None, None]
    if passed:
        tail = (1 - confidence_level) / 2
        bounds = np.quantile(replicates[valid], [tail, 1 - tail]).tolist()
    if not passed:
        status = "bootstrap_gate_failed"
    else:
        status = "ok" if n_valid == len(replicates) else "conditional_approximate"
    return {
        "estimate": float(estimate),
        "ci_lower": bounds[0],
        "ci_upper": bounds[1],
        "n_bootstrap_success": n_valid,
        "n_bootstrap_failed": len(replicates) - n_valid,
        "bootstrap_success_fraction": fraction,
        "conditional_on_defined_replicates": n_valid < len(replicates),
        "status": status,
    }


def _summarize(
    points: dict[str, np.ndarray],
    replicates: dict[str, np.ndarray],
    times: np.ndarray,
    confidence_level: float,
    minimum_success_fraction: float,
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for name in SCALAR_METRICS:
        result[name] = _interval_summary(
            float(points[name]), replicates[name], confidence_level, minimum_success_fraction
        )
    for name in CURVE_METRICS:
        result[name] = [
            {
                "time_months": float(time),
                **_interval_summary(
                    float(points[name][j]),
                    replicates[name][:, j],
                    confidence_level,
                    minimum_success_fraction,
                ),
            }
            for j, time in enumerate(times)
        ]
    summaries = [result[name] for name in SCALAR_METRICS]
    summaries.extend(item for name in CURVE_METRICS for item in result[name])
    result["reportable"] = all(item["status"] != "bootstrap_gate_failed" for item in summaries)
    if not result["reportable"]:
        result["status"] = "bootstrap_gate_failed"
    else:
        result["status"] = (
            "conditional_approximate"
            if any(item["status"] == "conditional_approximate" for item in summaries)
            else "ok"
        )
    return result


def evaluate_models(
    y_train: np.ndarray,
    y_test: np.ndarray,
    predictions: Mapping[str, Mapping[str, np.ndarray]],
    times: np.ndarray,
    config: Any,
    rng: np.random.Generator,
    *,
    bootstrap_strata: np.ndarray | None = None,
) -> dict[str, Any]:
    """Evaluate fixed fits with the same bootstrap indices for every model.

    Inputs contain no patient identifiers. The caller owns the only read of the
    frozen test split. This function never loads data, fits predictive models,
    or changes the grid in response to test outcomes. Intervals condition on the
    fitted models and the training censoring estimate; they are not refit CIs.
    Supply boolean PRIMARY-endpoint event strata from the frozen split even for
    an OS sensitivity. This keeps the original stratum counts in each replicate.
    """
    validate_survival_labels(y_train, len(y_train))
    try:
        validate_survival_labels(y_test, len(y_test))
    except ValueError as error:
        raise EvaluationSupportError(
            "Held-out survival labels are invalid or contain no observed events",
            {"n_test": len(y_test)},
        ) from error
    settings = _settings(config)
    times = np.asarray(times, dtype=float)
    expected_times = training_time_grid(y_train, settings)
    if not np.array_equal(times, expected_times):
        raise ValueError("Evaluation times must equal the prespecified training-only time grid")
    n_replicates = config_value(settings, "bootstrap_replicates")
    confidence_level = float(config_value(settings, "confidence_level"))
    minimum_fraction = float(config_value(settings, "minimum_bootstrap_success_fraction"))
    primary_metric = config_value(settings, "primary_metric")
    if (
        not isinstance(n_replicates, int)
        or isinstance(n_replicates, bool)
        or n_replicates < 2
        or not 0 < confidence_level < 1
        or not 0 < minimum_fraction <= 1
        or primary_metric not in SCALAR_METRICS
    ):
        raise ValueError("Invalid bootstrap settings or unknown primary scalar metric")
    if not predictions:
        raise ValueError("At least one set of model predictions is required")
    if not all(isinstance(name, str) for name in predictions):
        raise ValueError("Model identifiers must be strings")
    stratum_indices = None
    stratum_counts = None
    if bootstrap_strata is not None:
        bootstrap_strata = np.asarray(bootstrap_strata)
        if bootstrap_strata.shape != (len(y_test),) or bootstrap_strata.dtype.kind != "b":
            raise ValueError("bootstrap_strata must be aligned boolean primary-event labels")
        levels = np.unique(bootstrap_strata)
        stratum_indices = [np.flatnonzero(bootstrap_strata == level) for level in levels]
        stratum_counts = {
            str(bool(level)): len(indices)
            for level, indices in zip(levels, stratum_indices, strict=True)
        }

    event_name, time_name = y_test.dtype.names
    event, time = y_test[event_name], y_test[time_name]
    diagnostics = {
        "n_test": len(y_test),
        "n_test_events": int(event.sum()),
        "n_events_by_first_grid_time": int((event & (time <= times[0])).sum()),
        "n_at_risk_beyond_last_grid_time": int((time > times[-1]).sum()),
        "time_grid_months": times.tolist(),
    }
    if (
        time.min() >= times[0]
        or time.max() <= times[-1]
        or diagnostics["n_events_by_first_grid_time"] == 0
    ):
        raise EvaluationSupportError(
            "Held-out data do not support the prespecified grid; the grid was not changed",
            diagnostics,
        )
    eval_y, n_capped = _capped_test_labels(y_test, times[-1])
    train_censoring = CensoringDistributionEstimator().fit(y_train)
    supported_times = np.concatenate((times, eval_y[time_name]))
    try:
        censor_survival = train_censoring.predict_proba(supported_times)
    except ValueError as error:
        raise EvaluationSupportError(
            "Evaluation requires censoring probabilities beyond training support", diagnostics
        ) from error
    if (censor_survival <= 0).any():
        raise EvaluationSupportError(
            "Training censoring survival is zero on evaluation support", diagnostics
        )

    validated: dict[str, dict[str, np.ndarray]] = {}
    for name in sorted(predictions):
        risk = np.asarray(predictions[name]["risk"], dtype=float)
        survival = np.asarray(predictions[name]["survival"], dtype=float)
        if risk.shape != (len(y_test),) or survival.shape != (len(y_test), len(times)):
            raise ValueError("Prediction shapes must match held-out rows and the shared time grid")
        if not np.isfinite(risk).all() or not np.isfinite(survival).all():
            raise ValueError("Model predictions must be finite")
        if ((survival < 0) | (survival > 1)).any() or (np.diff(survival, axis=1) > 0).any():
            raise ValueError(
                "Survival probabilities must be in [0, 1] and non-increasing over time"
            )
        validated[name] = {"risk": risk, "survival": survival}

    points = {
        name: _metric_values(y_train, eval_y, arrays["risk"], arrays["survival"], times)
        for name, arrays in validated.items()
    }
    failed_points = {
        name: [metric for metric, value in values.items() if not np.isfinite(value).all()]
        for name, values in points.items()
    }
    failed_points = {name: failed for name, failed in failed_points.items() if failed}
    if failed_points:
        raise EvaluationSupportError(
            "Full-sample metrics are undefined on the prespecified grid",
            {**diagnostics, "undefined_metrics": failed_points},
        )
    bootstraps = {
        name: {
            **{metric: np.full(n_replicates, np.nan) for metric in SCALAR_METRICS},
            **{metric: np.full((n_replicates, len(times)), np.nan) for metric in CURVE_METRICS},
        }
        for name in validated
    }
    for b in range(n_replicates):
        if stratum_indices is None:
            indices = rng.integers(0, len(eval_y), size=len(eval_y))
        else:
            indices = np.concatenate(
                [group[rng.integers(0, len(group), size=len(group))] for group in stratum_indices]
            )
        boot_y = eval_y[indices]
        for name, arrays in validated.items():
            values = _metric_values(
                y_train, boot_y, arrays["risk"][indices], arrays["survival"][indices], times
            )
            for metric, value in values.items():
                bootstraps[name][metric][b] = value

    models = {
        name: _summarize(points[name], bootstraps[name], times, confidence_level, minimum_fraction)
        for name in validated
    }
    comparisons = {}
    for clinical_name in validated:
        if not clinical_name.endswith("__clinical"):
            continue
        model_name = clinical_name.removesuffix("__clinical")
        molecular_name = f"{model_name}__clinical_molecular"
        if molecular_name not in validated:
            continue
        difference = {
            metric: points[molecular_name][metric] - points[clinical_name][metric]
            for metric in points[clinical_name]
        }
        paired_replicates = {
            metric: bootstraps[molecular_name][metric] - bootstraps[clinical_name][metric]
            for metric in points[clinical_name]
        }
        comparisons[model_name] = {
            "difference": "clinical_molecular_minus_clinical",
            "positive_favors_molecular": {
                "uno_c": True,
                "ibs": False,
                "mean_auc": True,
                "brier": False,
                "auc": True,
            },
            **_summarize(difference, paired_replicates, times, confidence_level, minimum_fraction),
        }
    reportable = all(item["reportable"] for item in [*models.values(), *comparisons.values()])
    status = "ok"
    if any(
        item["status"] == "conditional_approximate"
        for item in [*models.values(), *comparisons.values()]
    ):
        status = "conditional_approximate"
    if not reportable:
        status = "bootstrap_gate_failed"
    return {
        "schema_version": 1,
        "status": status,
        "reportable": reportable,
        "n_training": len(y_train),
        **diagnostics,
        "evaluation_horizon_months": float(times[-1]),
        "test_label_administrative_censoring_months": float(np.nextafter(times[-1], np.inf)),
        "n_test_followup_capped": n_capped,
        "censoring_estimator": "training_reverse_kaplan_meier",
        "censoring_assumption": "C independent of (T,X); marginal censoring weights",
        "minimum_censoring_survival_on_evaluation_support": float(censor_survival.min()),
        "auc_risk_definition": "1_minus_survival_at_each_grid_time",
        "primary_metric": primary_metric,
        "bootstrap": {
            "replicates": n_replicates,
            "confidence_level": confidence_level,
            "minimum_success_fraction": minimum_fraction,
            "method": "paired_test_patient_percentile",
            "sampling_design": "fixed_primary_event_strata"
            if stratum_indices is not None
            else "iid_approximation_for_stratified_split",
            "stratum_counts": stratum_counts,
            "models_refitted": False,
            "training_censoring_weights_fixed": True,
            "multiple_comparisons_adjusted": False,
        },
        "models": models,
        "paired_feature_comparisons": comparisons,
    }
