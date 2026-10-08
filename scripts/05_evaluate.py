"""Stage 05: freeze the analysis, read test once, and evaluate fixed predictions.

The local ignored test cache supports interruption recovery without reopening raw
partitions. Once plan.json exists, stages03/04 prohibit refitting or recalibration.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
from pathlib import Path

import numpy as np
import yaml

from brca.calibration import load_calibration, model_quantiles
from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT
from brca.conformal.qin import bootstrap_cox_pivots, predict_qin_intervals
from brca.conformal.survival import (
    ConformalCalibration,
    PredictionIntervals,
    coverage_summary,
    paired_interval_comparison,
    subgroup_coverage,
)
from brca.data.eda import outcome_labels
from brca.data.partitions import load_partition
from brca.data.preprocessing import TrainingPreprocessor
from brca.evaluation.metrics import EvaluationSupportError, evaluate_models
from brca.fitting import load_fit, load_fit_index
from brca.manifest import sha256_file, write_json_exclusive
from brca.pipeline import StageRun, named_rng, write_json


def _save_cache(path, payload):
    temporary = path.with_suffix(".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def frozen_test_cache(run, index, calibration):
    """The sole downstream reader of test; the plan is written BEFORE the read."""
    directory = run.outputs / "evaluation"
    directory.mkdir(parents=True, exist_ok=True)
    plan_path, cache_path = directory / "plan.json", directory / "test_cache.pkl"
    protocol = {
        "config_sha256": run.metadata["config_sha256"],
        "split_sha256": index["split_sha256"],
        "fit_index_sha256": calibration["fit_index_sha256"],
        "calibration_sha256": sha256_file(run.outputs / "models/conformal_calibration.json"),
        "analysis_source_sha256": {
            str(path.relative_to(PROJECT_ROOT)): sha256_file(path)
            for path in sorted(
                [
                    PROJECT_ROOT / "scripts/05_evaluate.py",
                    *[
                        p
                        for p in (PROJECT_ROOT / "src/brca").rglob("*.py")
                        if "viz" not in p.parts and p.name != "reporting.py"
                    ],
                ]
            )
        },
    }
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
        if any(plan[key] != value for key, value in protocol.items()):
            raise ValueError(
                "The opened test analysis cannot change models, calibration or config."
            )
        if not cache_path.exists() or not plan.get("test_cache_sha256"):
            raise ValueError("Test opening was interrupted before a verified cache was saved.")
        if sha256_file(cache_path) != plan["test_cache_sha256"]:
            raise ValueError("Frozen test cache changed.")
        return pickle.loads(cache_path.read_bytes())
    plan = {
        **protocol,
        "source_sha256": run.metadata["source_sha256"],
        "source_archive": run.metadata["source_archive"],
        "test_partition_sha256": json.loads(
            (run.processed / "partitions/checksums.json").read_text()
        )["partitions_sha256"]["test"],
        "primary_endpoint": run.config.outcome.primary_endpoint,
        "primary_metric": run.config.evaluation.primary_metric,
        "time_grids": index["time_grids"],
        "status": "opened",
        "slurm_job_id": run.metadata["slurm_job_id"],
    }
    write_json_exclusive(plan_path, plan)
    test = load_partition(run.config, "test", allow_test=True).drop(columns="patient_id")
    _save_cache(cache_path, test)
    plan["test_cache_sha256"] = sha256_file(cache_path)
    plan["status"] = "cached"
    write_json(plan_path, plan)
    return test


def suppress_small_groups(groups, minimum):
    return {
        column: {
            level: (
                summary
                if summary["n_test"] >= minimum
                else {
                    "n_test": summary["n_test"],
                    "status": "suppressed_small_group",
                }
            )
            for level, summary in levels.items()
        }
        for column, levels in groups.items()
    }


def prediction_cache_binding(run, record):
    plan = json.loads((run.outputs / "evaluation/plan.json").read_text())
    return {
        "plan_sha256": sha256_file(run.outputs / "evaluation/plan.json"),
        "model_artifacts_sha256": record["artifacts_sha256"],
        "test_cache_sha256": plan["test_cache_sha256"],
    }


def run_qin_exploratory(run, index, qin_config_path):
    """Post hoc comparator; reuse the already-opened, hash-verified test cache."""
    cfg = run.config
    qin_config_path = Path(qin_config_path).resolve()
    qin_config = yaml.safe_load(qin_config_path.read_text())
    original_plan_path = run.outputs / "evaluation/plan.json"
    original_plan = json.loads(original_plan_path.read_text())
    fit_index_path = run.outputs / "models/fit_index.json"
    if (
        sha256_file(fit_index_path) != original_plan["fit_index_sha256"]
        or sha256_file(run.processed / "partitions/train.csv") != index["training_sha256"]
        or run.metadata["config_sha256"] != original_plan["config_sha256"]
    ):
        raise ValueError("Qin comparator inputs differ from the frozen original analysis")
    cache_path = run.outputs / "evaluation/test_cache.pkl"
    if sha256_file(cache_path) != original_plan["test_cache_sha256"]:
        raise ValueError("Frozen original test cache changed")
    directory = run.outputs / "qin_exploratory"
    directory.mkdir(parents=True, exist_ok=True)
    plan_path = directory / "plan.json"
    protocol = {
        "method": "qin_2025_bootstrap_cox_two_sided",
        "status": "post_hoc_exploratory_test_reuse",
        "config_sha256": run.metadata["config_sha256"],
        "training_sha256": index["training_sha256"],
        "original_plan_sha256": sha256_file(original_plan_path),
        "original_test_cache_sha256": sha256_file(cache_path),
        "fit_index_sha256": sha256_file(fit_index_path),
        "qin_config_sha256": sha256_file(qin_config_path),
        "analysis_source_sha256": {
            str(path.relative_to(PROJECT_ROOT)): sha256_file(path)
            for path in sorted(
                [
                    PROJECT_ROOT / "scripts/05_evaluate.py",
                    *[
                        p
                        for p in (PROJECT_ROOT / "src/brca").rglob("*.py")
                        if "viz" not in p.parts and p.name != "reporting.py"
                    ],
                ]
            )
        },
        "alphas": cfg.conformal.alphas,
        "horizons_months": cfg.conformal.horizons_months,
    }
    if plan_path.exists():
        if json.loads(plan_path.read_text()) != protocol:
            raise ValueError("Qin exploratory plan changed after opening test")
    else:
        write_json_exclusive(plan_path, protocol)
    # Only stage 05 accesses the existing test cache; no model choice uses its labels.
    test = pickle.loads(cache_path.read_bytes())
    train = load_partition(cfg, "train").drop(columns="patient_id")
    results = {}
    for endpoint in [cfg.outcome.primary_endpoint, *cfg.outcome.sensitivity_endpoints]:
        y_train, y_test = outcome_labels(train, endpoint), outcome_labels(test, endpoint)
        for feature_set in ("clinical", "clinical_molecular"):
            model_key = f"{endpoint}__cox__{feature_set}"
            model = load_fit(index["models"][model_key])
            preprocessor = TrainingPreprocessor.load(Path(index["preprocessors"][feature_set]))
            X_train, X_test = preprocessor.transform(train), preprocessor.transform(test)
            bootstrap = bootstrap_cox_pivots(
                X_train,
                y_train,
                model,
                cfg.cox,
                replicates=qin_config["bootstrap_replicates"],
                minimum_success_fraction=qin_config["minimum_bootstrap_success_fraction"],
                minimum_censoring_survival=qin_config["minimum_censoring_survival_for_event_pool"],
                rng=named_rng(cfg.seed, "qin_bootstrap", endpoint, feature_set),
            )
            for alpha in cfg.conformal.alphas:
                lower_pivot, upper_pivot = bootstrap.thresholds(alpha)
                for horizon in cfg.conformal.horizons_months:
                    key = f"{model_key}__alpha{alpha:g}__horizon{horizon:g}__qin_bootstrap_cox"
                    intervals = predict_qin_intervals(
                        model, X_test, bootstrap, alpha=alpha, horizon=horizon
                    )
                    coverage_kwargs = {
                        "bootstrap_replicates": cfg.evaluation.bootstrap_replicates,
                        "seed": int(
                            named_rng(cfg.seed, "coverage_bootstrap", endpoint).integers(
                                0, np.iinfo(np.uint32).max
                            )
                        ),
                        "confidence_level": cfg.evaluation.confidence_level,
                        "min_censoring_survival": cfg.conformal.min_censoring_survival,
                        "minimum_bootstrap_success_fraction": (
                            cfg.evaluation.minimum_bootstrap_success_fraction
                        ),
                        "bootstrap_strata": test["event"].to_numpy(bool),
                    }
                    groups = subgroup_coverage(
                        y_test,
                        intervals,
                        y_train,
                        groups={
                            column: test[column].to_numpy()
                            for column in cfg.conformal.coverage_subgroups
                        },
                        expected_levels={
                            column: sorted(train[column].dropna().unique().tolist())
                            for column in cfg.conformal.coverage_subgroups
                        },
                        **coverage_kwargs,
                    )
                    results[key] = {
                        "endpoint": endpoint,
                        "feature_set": feature_set,
                        "model": "cox",
                        "method": "qin_bootstrap_cox",
                        "alpha": alpha,
                        "horizon_months": horizon,
                        "bootstrap_pivot_quantiles": [lower_pivot, upper_pivot],
                        "bootstrap_diagnostics": bootstrap.diagnostics,
                        "marginal": coverage_summary(y_test, intervals, y_train, **coverage_kwargs),
                        "subgroups": suppress_small_groups(groups, cfg.eda.minimum_group_size),
                    }
            print(f"Qin bootstrap evaluated: {model_key}", flush=True)
    path = directory / "coverage.json"
    write_json(path, results)
    run.finish([plan_path, path], status="post_hoc_exploratory_complete", n_test=len(test))
    print(f"Qin aggregate report: {path}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--qin-exploratory", action="store_true")
    modes.add_argument("--reliability-exploratory", action="store_true")
    parser.add_argument("--qin-config", default=str(PROJECT_ROOT / "config/qin_exploratory.yaml"))
    parser.add_argument(
        "--reliability-config", default=str(PROJECT_ROOT / "config/uncertainty_reliability.yaml")
    )
    args = parser.parse_args()
    require_compute_node("stage 05 (evaluate)")
    run = StageRun(
        "05_uncertainty_reliability"
        if args.reliability_exploratory
        else ("05_qin_exploratory" if args.qin_exploratory else "05_evaluate"),
        args.config,
        allow_dirty=args.allow_dirty,
    )
    cfg = run.config
    index = load_fit_index(run)
    if args.reliability_exploratory:
        from brca.reliability import (
            evaluate_reliability,
            freeze_reliability_plan,
            prepare_reliability,
        )

        settings = freeze_reliability_plan(run, index, args.reliability_config)
        prepared = prepare_reliability(run, index, settings)
        cache = run.outputs / "evaluation/test_cache.pkl"
        protocol = json.loads((run.outputs / "uncertainty_reliability/plan.json").read_text())
        if sha256_file(cache) != protocol["test_cache_sha256"]:
            raise ValueError("Frozen original test cache changed during preparation")
        # This is the sole reader; preparation and all decisions are complete first.
        test = pickle.loads(cache.read_bytes())
        evaluate_reliability(run, index, prepared, test)
        return
    if args.qin_exploratory:
        run_qin_exploratory(run, index, args.qin_config)
        return
    calibration = load_calibration(run, index)
    train = load_partition(cfg, "train")
    test = frozen_test_cache(run, index, calibration)
    predictions, base_quantiles = {}, {}
    for key, record in index["models"].items():
        path = run.outputs / "evaluation" / f"predictions__{key}.pkl"
        metadata_path = path.with_suffix(".json")
        binding = prediction_cache_binding(run, record)
        if path.exists():
            metadata = json.loads(metadata_path.read_text())
            if metadata["binding"] != binding or metadata["sha256"] != sha256_file(path):
                raise ValueError("Prediction cache changed or belongs to another analysis.")
            cached = pickle.loads(path.read_bytes())
        else:
            model = load_fit(record)
            preprocessor = TrainingPreprocessor.load(
                Path(index["preprocessors"][record["feature_set"]])
            )
            X = preprocessor.transform(test)
            times = np.asarray(index["time_grids"][record["endpoint"]])
            cached = {
                "risk": model.predict_risk(X),
                "survival": model.predict_survival_function(X, times),
                "quantiles": model_quantiles(model, X, cfg.conformal.alphas),
            }
            _save_cache(path, cached)
            write_json(metadata_path, {"binding": binding, "sha256": sha256_file(path)})
        predictions[key] = {name: cached[name] for name in ("risk", "survival")}
        base_quantiles[key] = cached["quantiles"]
    performance = {}
    for endpoint in [cfg.outcome.primary_endpoint, *cfg.outcome.sensitivity_endpoints]:
        endpoint_predictions = {
            key.removeprefix(f"{endpoint}__"): value
            for key, value in predictions.items()
            if index["models"][key]["endpoint"] == endpoint
        }
        try:
            performance[endpoint] = evaluate_models(
                outcome_labels(train, endpoint),
                outcome_labels(test, endpoint),
                endpoint_predictions,
                np.asarray(index["time_grids"][endpoint]),
                cfg,
                named_rng(cfg.seed, "performance_bootstrap", endpoint),
                bootstrap_strata=test["event"].to_numpy(bool),
            )
        except EvaluationSupportError as error:
            failure_path = run.outputs / "metrics/evaluation_support_failure.json"
            write_json(
                failure_path,
                {"endpoint": endpoint, "reason": str(error), "diagnostics": error.diagnostics},
            )
            run.finish([failure_path], status="unsupported_prespecified_grid")
            raise
        print(f"Performance evaluated: {endpoint}", flush=True)
    performance_path = run.outputs / "metrics/performance.json"
    write_json(performance_path, performance)
    coverage, interval_objects = {}, {}
    states = dict(calibration["states"])
    # Native Bayesian predictive intervals are evaluated alongside conformal intervals.
    for key, record in index["models"].items():
        if record["model"] != "weibull":
            continue
        for alpha in cfg.conformal.alphas:
            for horizon in cfg.conformal.horizons_months:
                states[f"{key}__alpha{alpha:g}__horizon{horizon:g}__model_predictive"] = {
                    "model_key": key,
                    "endpoint": record["endpoint"],
                    "feature_set": record["feature_set"],
                    "model": record["model"],
                    "registered_method": "model_predictive",
                    "calibration": {"alpha": alpha, "horizon": horizon},
                }
    for key, state in states.items():
        endpoint = state["endpoint"]
        parameters = state["calibration"]
        alpha, horizon = parameters["alpha"], parameters["horizon"]
        lower, upper = np.minimum(base_quantiles[state["model_key"]][alpha], horizon).T
        if state["registered_method"] == "model_predictive":
            intervals = PredictionIntervals(lower, upper, alpha, 0, horizon)
        else:
            intervals = ConformalCalibration.from_dict(parameters).predict(lower, upper)
        interval_objects[key] = intervals
        y_test, y_train = outcome_labels(test, endpoint), outcome_labels(train, endpoint)
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
        coverage[key] = {
            **{k: state[k] for k in ("endpoint", "model", "feature_set", "registered_method")},
            "marginal": coverage_summary(y_test, intervals, y_train, **kwargs),
            "subgroups": suppress_small_groups(groups, cfg.eda.minimum_group_size),
        }
    coverage_path = run.outputs / "metrics/conformal_coverage.json"
    write_json(coverage_path, coverage)
    paired_intervals = {}
    for key, state in states.items():
        if state["feature_set"] != "clinical":
            continue
        molecular = key.replace("__clinical__", "__clinical_molecular__")
        endpoint = state["endpoint"]
        paired_intervals[key] = paired_interval_comparison(
            outcome_labels(test, endpoint),
            interval_objects[key],
            interval_objects[molecular],
            outcome_labels(train, endpoint),
            bootstrap_replicates=cfg.evaluation.bootstrap_replicates,
            seed=int(
                named_rng(cfg.seed, "coverage_bootstrap", endpoint).integers(
                    0, np.iinfo(np.uint32).max
                )
            ),
            confidence_level=cfg.evaluation.confidence_level,
            min_censoring_survival=cfg.conformal.min_censoring_survival,
            minimum_bootstrap_success_fraction=cfg.evaluation.minimum_bootstrap_success_fraction,
            bootstrap_strata=test["event"].to_numpy(bool),
        )
    comparisons_path = run.outputs / "metrics/rq2_feature_set_comparison.json"
    write_json(
        comparisons_path,
        {
            "performance": {
                endpoint: value["paired_feature_comparisons"]
                for endpoint, value in performance.items()
            },
            "intervals": paired_intervals,
        },
    )
    artifacts = [
        performance_path,
        coverage_path,
        comparisons_path,
        run.outputs / "evaluation/plan.json",
        *sorted((run.outputs / "evaluation").glob("predictions__*.json")),
    ]
    reportable = all(report["reportable"] for report in performance.values())
    qualified = any(report["status"] != "ok" for report in performance.values())
    status = "complete" if not qualified else "complete_with_qualified_intervals"
    if not reportable:
        status = "bootstrap_gate_failed"
    run.finish(artifacts, status=status, n_test=len(test), reportable=reportable)
    print(f"Final aggregate metrics: {run.outputs / 'metrics'}")


if __name__ == "__main__":
    main()
