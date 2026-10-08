"""Isolated post hoc survival-uncertainty reliability experiment.

One-sided score envelopes use frozen fits and the original calibration partition.
A fixed-ridge Qin comparator diagnoses numerical robustness, not clinical validity.
Only the stage-05 caller reads the frozen test cache; this module receives its frame.
"""

from __future__ import annotations

import json
import warnings
from collections import Counter
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml
from scipy.linalg import LinAlgWarning
from sklearn.exceptions import ConvergenceWarning
from sksurv.linear_model import CoxPHSurvivalAnalysis

from brca.calibration import interval_key, load_calibration
from brca.config import PROJECT_ROOT
from brca.conformal.lower import LowerBoundCalibration, calibrate_lower_bound
from brca.conformal.qin import (
    QinBootstrap,
    event_pool_weights,
    inverse_survival_times,
    predict_qin_intervals,
)
from brca.conformal.qin_ridge import (
    QinRidgeBootstrapError,
    bootstrap_ridge_cox_pivots,
    fit_qin_ridge_model,
)
from brca.conformal.survival import (
    CensoringDistribution,
    ConformalCalibration,
    PredictionIntervals,
    coverage_summary,
    paired_interval_comparison,
    subgroup_coverage,
)
from brca.data.eda import outcome_labels
from brca.data.partitions import load_partition
from brca.data.preprocessing import TrainingPreprocessor
from brca.evaluation.metrics import evaluate_models
from brca.fitting import load_fit
from brca.manifest import sha256_file, write_json_exclusive
from brca.models.cox import CoxPHModel
from brca.pipeline import named_rng, write_json


def reliability_directory(run):
    return run.outputs / "uncertainty_reliability"


def freeze_reliability_plan(run, index, settings_path):
    """Bind this separate experiment to original inputs and current statistical code."""
    settings_path = Path(settings_path).resolve()
    settings = yaml.safe_load(settings_path.read_text())
    if (
        set(settings)
        != {
            "ridge_alpha",
            "ridge_bootstrap_replicates",
            "diagnostic_bootstrap_replicates",
            "report_alpha",
        }
        or not np.isfinite(settings["ridge_alpha"])
        or settings["ridge_alpha"] <= 0
        or any(
            type(settings[key]) is not int or settings[key] < 2
            for key in ("ridge_bootstrap_replicates", "diagnostic_bootstrap_replicates")
        )
        or settings["report_alpha"] not in run.config.conformal.alphas
    ):
        raise ValueError("Invalid reliability settings")
    original_path = run.outputs / "evaluation/plan.json"
    original = json.loads(original_path.read_text())
    index_path = run.outputs / "models/fit_index.json"
    cache = run.outputs / "evaluation/test_cache.pkl"
    calibration = load_calibration(run, index)
    if (
        sha256_file(index_path) != original["fit_index_sha256"]
        or run.metadata["config_sha256"] != original["config_sha256"]
        or sha256_file(run.processed / "partitions/train.csv") != index["training_sha256"]
        or sha256_file(cache) != original["test_cache_sha256"]
        or sha256_file(run.outputs / "models/conformal_calibration.json")
        != original["calibration_sha256"]
        or sha256_file(run.processed / "partitions/calibration.csv")
        != calibration["calibration_partition_sha256"]
    ):
        raise ValueError("Reliability experiment inputs differ from the original frozen analysis")
    sources = sorted(
        [
            PROJECT_ROOT / "scripts/05_evaluate.py",
            *[
                path
                for path in (PROJECT_ROOT / "src/brca").rglob("*.py")
                if "viz" not in path.parts and path.name != "reporting.py"
            ],
        ]
    )
    protocol = {
        "status": "post_hoc_exploratory_test_reuse",
        "methods": ["conservative_lower_bound", "qin_bootstrap_ridge_cox"],
        "config_sha256": run.metadata["config_sha256"],
        "settings_sha256": sha256_file(settings_path),
        "settings": settings,
        "lock_sha256": sha256_file(PROJECT_ROOT / "uv.lock"),
        "fit_index_sha256": sha256_file(index_path),
        "original_plan_sha256": sha256_file(original_path),
        "test_cache_sha256": original["test_cache_sha256"],
        "split_sha256": index["split_sha256"],
        "partition_sha256": {
            name: sha256_file(run.processed / f"partitions/{name}.csv")
            for name in ("train", "calibration")
        },
        "reference_artifacts_sha256": {
            str(path): sha256_file(path)
            for path in (
                run.outputs / "models/conformal_calibration.json",
                run.outputs / "metrics/conformal_coverage.json",
                run.outputs / "metrics/performance.json",
                run.outputs / "qin_exploratory/coverage.json",
                run.outputs / "qin_exploratory/plan.json",
            )
        },
        "analysis_source_sha256": {
            str(path.relative_to(PROJECT_ROOT)): sha256_file(path) for path in sources
        },
        "alphas": run.config.conformal.alphas,
        "horizons_months": run.config.conformal.horizons_months,
    }
    directory = reliability_directory(run)
    plan_path = directory / "plan.json"
    directory.mkdir(parents=True, exist_ok=True)
    if plan_path.exists():
        if json.loads(plan_path.read_text()) != protocol:
            raise ValueError("Reliability protocol changed; preserve the opened experiment")
    else:
        write_json_exclusive(plan_path, protocol)
    return settings


def _ess(weights):
    return float(weights.sum() ** 2 / np.square(weights).sum()) if len(weights) else 0.0


def training_support(y_train, horizons, minimum_survival):
    """Aggregate support audit: never pick a horizon by test results."""
    event, time = y_train["event"], y_train["time"]
    censoring = CensoringDistribution.fit(y_train)
    _, _, full = event_pool_weights(y_train, minimum_censoring_survival=0.0)
    by_horizon = {}
    for horizon in horizons:
        g = float(censoring.survival_ge(np.array([horizon]))[0])
        before = event & (time < horizon)
        through = time >= horizon
        event_weights = 1 / censoring.survival_ge(time[before])
        supported = bool(np.isfinite(g) and g > 0)
        known_weights = (
            np.concatenate([event_weights, np.full(through.sum(), 1 / g)]) if supported else None
        )
        by_horizon[str(horizon)] = {
            "horizon_months": horizon,
            "censoring_survival_at_horizon": g if np.isfinite(g) else None,
            "ipcw_support_status": "supported"
            if supported and g >= minimum_survival
            else "unsupported",
            "n_events_before_horizon": int(before.sum()),
            "n_observed_through_horizon": int(through.sum()),
            "pre_horizon_event_effective_sample_size": _ess(event_weights),
            "restricted_known_target_effective_sample_size": _ess(known_weights)
            if supported
            else None,
            "weighted_horizon_mass": float((through.sum() / g) / known_weights.sum())
            if supported and known_weights.sum() > 0
            else None,
            "interpretation": "diagnostic_only_not_an_alternative_calibration",
        }
    return {"n_train": len(y_train), "full_event_pool": full, "horizons": by_horizon}


def diagnose_unpenalized_bootstrap(X, y, original_model, config, *, replicates, rng):
    """Classify the old numerical failures using only training resamples."""
    active = X.loc[:, list(original_model.active_feature_names_)]
    indices, probabilities, _ = event_pool_weights(y, minimum_censoring_survival=0.0)
    counters, reasons = Counter(), Counter()
    for _ in range(replicates):
        sampled = rng.integers(0, len(X), size=len(X))
        rng.choice(indices, p=probabilities)  # Preserve the original Qin RNG sequence.
        values = active.iloc[sampled].to_numpy()
        constant = bool(np.any(np.ptp(values, axis=0) == 0))
        deficient = bool(np.linalg.matrix_rank(values - values.mean(axis=0)) < values.shape[1])
        counters["constant_column_resamples"] += int(constant)
        counters["rank_deficient_resamples"] += int(deficient)
        try:
            with warnings.catch_warnings():
                for category in (ConvergenceWarning, LinAlgWarning, RuntimeWarning):
                    warnings.simplefilter("error", category)
                model = CoxPHSurvivalAnalysis(
                    alpha=config.alpha, ties=config.ties, n_iter=config.n_iter, tol=config.tol
                ).fit(active.iloc[sampled], y[sampled])
            if not np.isfinite(model.coef_).all():
                raise ValueError("Nonfinite coefficients")
            counters["successful"] += 1
        except (
            ArithmeticError,
            ConvergenceWarning,
            LinAlgWarning,
            RuntimeWarning,
            ValueError,
            RuntimeError,
            np.linalg.LinAlgError,
        ) as error:
            counters["failed"] += 1
            counters["failed_with_constant_column"] += int(constant)
            counters["failed_with_rank_deficiency"] += int(deficient)
            reasons[type(error).__name__] += 1
    return {"requested": replicates, "counts": dict(counters), "failure_reasons": dict(reasons)}


def prepare_reliability(run, index, settings):
    """Only train and calibration are read here; no test-frame parameter exists."""
    cfg, directory = run.config, reliability_directory(run)
    state_path = directory / "preparation.json"
    if state_path.exists():
        state = json.loads(state_path.read_text())
        binding = json.loads((directory / "preparation_binding.json").read_text())
        if (
            binding["plan_sha256"] != sha256_file(directory / "plan.json")
            or binding["preparation_sha256"] != sha256_file(state_path)
            or any(
                sha256_file(path) != digest for path, digest in state["artifacts_sha256"].items()
            )
        ):
            raise ValueError("Reliability preparation cache changed")
        return state
    train = load_partition(cfg, "train")
    calibration = load_partition(cfg, "calibration")
    states, ridge, diagnostics, artifacts = {}, {}, {}, {}
    for endpoint in [cfg.outcome.primary_endpoint, *cfg.outcome.sensitivity_endpoints]:
        diagnostics[endpoint] = training_support(
            outcome_labels(train, endpoint),
            cfg.conformal.horizons_months,
            cfg.conformal.min_censoring_survival,
        )
    for key, record in index["models"].items():
        model = load_fit(record)
        pre = TrainingPreprocessor.load(Path(index["preprocessors"][record["feature_set"]]))
        X_cal = pre.transform(calibration)
        quantiles = model.predict_quantiles(X_cal, np.asarray(cfg.conformal.alphas))
        y_cal = outcome_labels(calibration, record["endpoint"])
        for column, alpha in enumerate(cfg.conformal.alphas):
            for horizon in cfg.conformal.horizons_months:
                fitted = calibrate_lower_bound(
                    y_cal,
                    quantiles[:, column],
                    alpha,
                    horizon=horizon,
                    strata=calibration["event"].to_numpy(bool),
                    expected_strata=(False, True),
                )
                states[interval_key(key, alpha, horizon, "conservative_lower_bound")] = {
                    "model_key": key,
                    "model": record["model"],
                    "endpoint": record["endpoint"],
                    "feature_set": record["feature_set"],
                    "calibration": fitted.to_dict(),
                }
        if record["model"] == "cox":
            X_train, y_train = pre.transform(train), outcome_labels(train, record["endpoint"])
            diagnostic = diagnose_unpenalized_bootstrap(
                X_train,
                y_train,
                model,
                cfg.cox,
                replicates=settings["diagnostic_bootstrap_replicates"],
                rng=named_rng(cfg.seed, "qin_bootstrap", record["endpoint"], record["feature_set"]),
            )
            ridge_config = SimpleNamespace(**{**vars(cfg.cox), "alpha": settings["ridge_alpha"]})
            ridge_model = fit_qin_ridge_model(
                X_train, y_train, ridge_config, rng=named_rng(cfg.seed, "qin_ridge_fit", key)
            )
            path = directory / f"ridge_model__{key}.pkl"
            ridge_model.save(path)
            artifacts[str(path)] = sha256_file(path)
            try:
                bootstrap = bootstrap_ridge_cox_pivots(
                    X_train,
                    y_train,
                    ridge_model,
                    ridge_config,
                    replicates=settings["ridge_bootstrap_replicates"],
                    minimum_censoring_survival=0.0,
                    rng=named_rng(
                        cfg.seed, "qin_bootstrap", record["endpoint"], record["feature_set"]
                    ),
                )
            except QinRidgeBootstrapError as error:
                failure = directory / f"failed_bootstrap__{key}.json"
                write_json(failure, error.diagnostics)
                run.finish([directory / "plan.json", failure], status="ridge_bootstrap_gate_failed")
                raise
            ridge[key] = {
                "path": str(path),
                "endpoint": record["endpoint"],
                "feature_set": record["feature_set"],
                "pivots": bootstrap.pivots.tolist(),
                "diagnostics": bootstrap.diagnostics,
                "unpenalized_pilot": diagnostic,
            }
        print(f"Reliability train/calibration prepared: {key}", flush=True)
    state = {
        "lower_states": states,
        "ridge": ridge,
        "training_support": diagnostics,
        "artifacts_sha256": artifacts,
        "n_train": len(train),
        "n_cal": len(calibration),
    }
    write_json(state_path, state)
    write_json(
        directory / "preparation_binding.json",
        {
            "plan_sha256": sha256_file(directory / "plan.json"),
            "preparation_sha256": sha256_file(state_path),
        },
    )
    return state


def _summary(run, train, test, endpoint, intervals):
    cfg = run.config
    kwargs = {
        "bootstrap_replicates": cfg.evaluation.bootstrap_replicates,
        "seed": int(
            named_rng(cfg.seed, "coverage_bootstrap", endpoint).integers(0, np.iinfo(np.uint32).max)
        ),
        "confidence_level": cfg.evaluation.confidence_level,
        "min_censoring_survival": cfg.conformal.min_censoring_survival,
        "minimum_bootstrap_success_fraction": cfg.evaluation.minimum_bootstrap_success_fraction,
        "bootstrap_strata": test["event"].to_numpy(bool),
    }
    y_train, y_test = outcome_labels(train, endpoint), outcome_labels(test, endpoint)
    groups = subgroup_coverage(
        y_test,
        intervals,
        y_train,
        groups={column: test[column].to_numpy() for column in cfg.conformal.coverage_subgroups},
        expected_levels={
            column: sorted(train[column].dropna().unique().tolist())
            for column in cfg.conformal.coverage_subgroups
        },
        **kwargs,
    )
    suppressed = {
        column: {
            level: value
            if value["n_test"] >= cfg.eda.minimum_group_size
            else {"n_test": value["n_test"], "status": "suppressed_small_group"}
            for level, value in levels.items()
        }
        for column, levels in groups.items()
    }
    return {
        "marginal": coverage_summary(y_test, intervals, y_train, **kwargs),
        "subgroups": suppressed,
    }, kwargs


def _paired_methods(y_test, reference, candidate, y_train, kwargs, *, direction):
    comparison = paired_interval_comparison(y_test, reference, candidate, y_train, **kwargs)
    comparison["direction"] = direction
    comparison["width_interpretation"] = "positive_difference_means_wider_candidate_intervals"
    comparison["n_cal_reference"] = comparison.pop("n_cal_clinical")
    comparison["n_cal_candidate"] = comparison.pop("n_cal_molecular")
    return comparison


def _verify_reference(y_test, intervals, y_train, saved):
    """Reject comparisons if reconstructed historical point statistics differ."""
    actual = coverage_summary(
        y_test,
        intervals,
        y_train,
        bootstrap_replicates=0,
        seed=0,
        confidence_level=0.95,
        min_censoring_survival=0.05,
        minimum_bootstrap_success_fraction=0.9,
    )
    for metric in (
        "observable_coverage_lower",
        "observable_coverage_upper",
        "width_median_months",
        "width_q25_months",
        "width_q75_months",
        "full_support_fraction",
    ):
        if not np.isclose(actual[metric], saved[metric], atol=1e-9, rtol=1e-10):
            raise ValueError(f"Historical interval reconstruction changed: {metric}")


def _contains_bootstrap_failure(value):
    if isinstance(value, dict):
        return value.get("status") == "bootstrap_gate_failed" or any(
            _contains_bootstrap_failure(item) for item in value.values()
        )
    if isinstance(value, list):
        return any(_contains_bootstrap_failure(item) for item in value)
    return False


def evaluate_reliability(run, index, state, test):
    """Evaluate all fixed candidates once; retain paired differences and all subgroups."""
    cfg, directory = run.config, reliability_directory(run)
    train = load_partition(cfg, "train")
    original_calibration = load_calibration(run, index)["states"]
    original_qin = json.loads((run.outputs / "qin_exploratory/coverage.json").read_text())
    original_coverage = json.loads((run.outputs / "metrics/conformal_coverage.json").read_text())
    reports, comparisons, feature_pairs, interval_objects, ridge_predictions = {}, {}, {}, {}, {}
    probabilities = sorted(
        {p for alpha in cfg.conformal.alphas for p in (alpha, alpha / 2, 1 - alpha / 2)}
    )
    for model_key, record in index["models"].items():
        endpoint, feature = record["endpoint"], record["feature_set"]
        model = load_fit(record)
        pre = TrainingPreprocessor.load(Path(index["preprocessors"][feature]))
        X = pre.transform(test)
        quantiles = model.predict_quantiles(X, np.asarray(probabilities))
        y_train, y_test = outcome_labels(train, endpoint), outcome_labels(test, endpoint)
        ridge_model = None
        if model_key in state["ridge"]:
            item = state["ridge"][model_key]
            ridge_model = CoxPHModel.load(item["path"])
            bootstrap = QinBootstrap(np.asarray(item["pivots"]), item["diagnostics"])
            times = np.asarray(index["time_grids"][endpoint])
            ridge_predictions.setdefault(endpoint, {})[f"cox_ridge__{feature}"] = {
                "risk": ridge_model.predict_risk(X),
                "survival": ridge_model.predict_survival_function(X, times),
            }
        for alpha in cfg.conformal.alphas:
            for horizon in cfg.conformal.horizons_months:
                key = interval_key(model_key, alpha, horizon, "conservative_lower_bound")
                fitted = LowerBoundCalibration.from_dict(state["lower_states"][key]["calibration"])
                lower = fitted.predict(quantiles[:, probabilities.index(alpha)])
                interval_objects[key] = lower
                reference_key = interval_key(model_key, alpha, horizon, "conservative")
                reference_state = original_calibration[reference_key]["calibration"]
                reference = ConformalCalibration.from_dict(reference_state).predict(
                    np.minimum(quantiles[:, probabilities.index(alpha / 2)], horizon),
                    np.minimum(quantiles[:, probabilities.index(1 - alpha / 2)], horizon),
                )
                _verify_reference(
                    y_test, reference, y_train, original_coverage[reference_key]["marginal"]
                )
                summary, kwargs = _summary(run, train, test, endpoint, lower)
                reports[key] = {
                    "endpoint": endpoint,
                    "feature_set": feature,
                    "model": record["model"],
                    "method": lower.method,
                    "calibration_sample_type": "held_out_patients",
                    **summary,
                }
                comparisons[key] = _paired_methods(
                    y_test,
                    reference,
                    lower,
                    y_train,
                    kwargs,
                    direction="lower_bound_minus_primary_conservative_same_model",
                )
                if ridge_model is not None:
                    new = replace(
                        predict_qin_intervals(
                            ridge_model, X, bootstrap, alpha=alpha, horizon=horizon
                        ),
                        method="qin_bootstrap_ridge_cox",
                    )
                    new_key = interval_key(model_key, alpha, horizon, new.method)
                    interval_objects[new_key] = new
                    summary, _ = _summary(run, train, test, endpoint, new)
                    reports[new_key] = {
                        "endpoint": endpoint,
                        "feature_set": feature,
                        "model": "cox_ridge",
                        "method": new.method,
                        "calibration_sample_type": "bootstrap_pivots_not_held_out_patients",
                        "bootstrap_diagnostics": bootstrap.diagnostics,
                        **summary,
                    }
                    old_key = interval_key(model_key, alpha, horizon, "qin_bootstrap_cox")
                    old = original_qin[old_key]
                    lo_pivot, hi_pivot = old["bootstrap_pivot_quantiles"]
                    survival = model.estimator_.predict_survival_function(
                        X.loc[:, model.active_feature_names_], return_array=True
                    )
                    knots = model.estimator_.unique_times_
                    old_intervals = PredictionIntervals(
                        np.minimum(inverse_survival_times(survival, knots, hi_pivot), horizon),
                        np.minimum(
                            inverse_survival_times(survival, knots, lo_pivot, strict=True), horizon
                        ),
                        alpha,
                        old["marginal"]["n_cal"],
                        horizon,
                        "qin_bootstrap_cox",
                    )
                    _verify_reference(y_test, old_intervals, y_train, old["marginal"])
                    comparisons[new_key] = _paired_methods(
                        y_test,
                        old_intervals,
                        new,
                        y_train,
                        kwargs,
                        direction="ridge_qin_minus_original_qin_cox",
                    )
        print(f"Reliability test evaluated: {model_key}", flush=True)
    for key, intervals in interval_objects.items():
        if "__clinical__" not in key:
            continue
        other = key.replace("__clinical__", "__clinical_molecular__")
        endpoint = reports[key]["endpoint"]
        # _summary uses the same deterministic test bootstrap seed for each endpoint.
        kwargs = {
            "bootstrap_replicates": cfg.evaluation.bootstrap_replicates,
            "seed": int(
                named_rng(cfg.seed, "coverage_bootstrap", endpoint).integers(
                    0, np.iinfo(np.uint32).max
                )
            ),
            "confidence_level": cfg.evaluation.confidence_level,
            "min_censoring_survival": cfg.conformal.min_censoring_survival,
            "minimum_bootstrap_success_fraction": cfg.evaluation.minimum_bootstrap_success_fraction,
            "bootstrap_strata": test["event"].to_numpy(bool),
        }
        feature_pairs[key] = paired_interval_comparison(
            outcome_labels(test, endpoint),
            intervals,
            interval_objects[other],
            outcome_labels(train, endpoint),
            **kwargs,
        )
    performance = {
        endpoint: evaluate_models(
            outcome_labels(train, endpoint),
            outcome_labels(test, endpoint),
            predictions,
            np.asarray(index["time_grids"][endpoint]),
            cfg,
            named_rng(cfg.seed, "performance_bootstrap", endpoint),
            bootstrap_strata=test["event"].to_numpy(bool),
        )
        for endpoint, predictions in ridge_predictions.items()
    }
    artifacts = []
    for name, payload in (
        ("coverage", reports),
        ("paired_method_comparisons", comparisons),
        ("paired_feature_comparisons", feature_pairs),
        ("ridge_performance", performance),
    ):
        path = directory / f"{name}.json"
        write_json(path, payload)
        artifacts.append(path)
    artifacts += [
        directory / "plan.json",
        directory / "preparation.json",
        directory / "preparation_binding.json",
        *map(Path, state["artifacts_sha256"]),
    ]
    failed = _contains_bootstrap_failure([reports, comparisons, feature_pairs, performance])
    run.finish(
        artifacts,
        status="bootstrap_gate_failed" if failed else "exploratory_complete",
        reportable=not failed,
        n_test=len(test),
        n_interval_cells=len(reports),
        historical_point_statistics_verified=True,
    )
    if failed:
        raise ValueError("Reliability evaluation bootstrap gate failed; inspect saved diagnostics")
    return directory
