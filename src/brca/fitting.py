"""Training-only orchestration and checked persistence of the registered model matrix."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from brca.config import PROJECT_ROOT
from brca.data.eda import outcome_labels
from brca.data.partitions import load_partition
from brca.data.preprocessing import TrainingPreprocessor
from brca.evaluation.metrics import training_time_grid
from brca.manifest import sha256_file
from brca.models.bayesian import (
    BayesianAFTModel,
    BayesianDiagnosticsError,
    BayesianInitializationError,
    PriorPredictiveError,
)
from brca.models.cox import CoxPHModel
from brca.models.rsf import RandomSurvivalForestModel
from brca.pipeline import StageRun, named_rng, write_json

FEATURE_SETS = ("clinical", "clinical_molecular")
MODEL_TYPES = {"cox": CoxPHModel, "rsf": RandomSurvivalForestModel, "weibull": BayesianAFTModel}


def verify_artifacts(record: dict) -> None:
    for path, expected in record["artifacts_sha256"].items():
        if not Path(path).is_file() or sha256_file(path) != expected:
            raise ValueError("A registered fit or preprocessing artifact changed.")


def load_fit(record: dict):
    verify_artifacts(record)
    return MODEL_TYPES[record["model"]].load(record["path"])


def load_fit_index(run: StageRun) -> dict:
    path = run.outputs / "models" / "fit_index.json"
    index = json.loads(path.read_text())
    if (
        index["config_sha256"] != run.metadata["config_sha256"]
        or index["split_sha256"] != sha256_file(run.processed / "splits.json")
        or not index.get("complete")
    ):
        raise ValueError("Fitted model matrix is incomplete or belongs to a different protocol.")
    verify_artifacts(index)
    for source, expected in index.get("fit_source_sha256", {}).items():
        if sha256_file(source) != expected:
            raise ValueError("Registered fitting or prediction source changed after fitting.")
    return index


def run_fitting(run: StageRun) -> Path:
    cfg = run.config
    if (run.outputs / "evaluation" / "plan.json").exists():
        raise ValueError("Test evaluation has opened; refitting this analysis is prohibited.")
    train = load_partition(cfg, "train")
    model_dir = run.outputs / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    index_path = model_dir / "fit_index.json"
    # Hash only code that determines fitted models; adding exports does not invalidate fits.
    fit_sources = [
        PROJECT_ROOT / "src/brca/fitting.py",
        PROJECT_ROOT / "src/brca/data/preprocessing.py",
        PROJECT_ROOT / "src/brca/data/eda.py",
        PROJECT_ROOT / "src/brca/pipeline.py",
    ]
    fit_sources += sorted((PROJECT_ROOT / "src/brca/models").glob("*.py"))
    protocol = {
        "config_sha256": run.metadata["config_sha256"],
        "split_sha256": sha256_file(run.processed / "splits.json"),
        "training_sha256": sha256_file(run.processed / "partitions/train.csv"),
        "fit_source_sha256": {str(p): sha256_file(p) for p in fit_sources},
    }
    index = {
        **protocol,
        "complete": False,
        "models": {},
        "sensitivities": {},
        "preprocessors": {},
        "time_grids": {},
        "artifacts_sha256": {},
    }
    if index_path.exists():
        index = json.loads(index_path.read_text())
        if any(index.get(key) != value for key, value in protocol.items()):
            raise ValueError(
                "Existing fits use different training inputs/config/code; do not overwrite."
            )
        verify_artifacts(index)
    matrices = {}
    for feature in FEATURE_SETS:
        path = model_dir / f"preprocessor__{feature}.pkl"
        if feature in index["preprocessors"]:
            preprocessor = TrainingPreprocessor.load(path)
        else:
            preprocessor = TrainingPreprocessor(cfg, feature).fit(train)
            preprocessor.save(path)
            index["preprocessors"][feature] = str(path)
            index["artifacts_sha256"][str(path)] = sha256_file(path)
        matrices[feature] = preprocessor.transform(train)
    tasks = []
    for endpoint in [cfg.outcome.primary_endpoint, *cfg.outcome.sensitivity_endpoints]:
        y = outcome_labels(train, endpoint)
        index["time_grids"][endpoint] = training_time_grid(y, cfg).tolist()
        for feature in FEATURE_SETS:
            for model_name in MODEL_TYPES:
                tasks.append((endpoint, feature, model_name, "default", "weibull", "models"))
    for feature in FEATURE_SETS:
        tasks.extend(
            [
                (
                    cfg.outcome.primary_endpoint,
                    feature,
                    "weibull",
                    "sensitivity",
                    "weibull",
                    "sensitivities",
                ),
                (
                    cfg.outcome.primary_endpoint,
                    feature,
                    "weibull",
                    "default",
                    "lognormal",
                    "sensitivities",
                ),
            ]
        )
    write_json(index_path, index)
    for endpoint, feature, name, prior, family, section in tasks:
        key = f"{endpoint}__{name}__{feature}"
        if section == "sensitivities":
            key += f"__{family}__{prior}"
        if key in index[section]:
            verify_artifacts(index[section][key])
            print(f"Reuse verified fit: {key}", flush=True)
            continue
        rng = named_rng(cfg.seed, "fit", endpoint, feature, name, prior, family)
        model = (
            BayesianAFTModel(cfg, rng, family=family, prior_variant=prior)
            if name == "weibull"
            else MODEL_TYPES[name](getattr(cfg, name), rng=rng)
        )
        path = model_dir / f"{key}{'.nc' if name == 'weibull' else '.pkl'}"
        metric_path = (
            run.outputs
            / "metrics"
            / f"{'bayes_diagnostics' if name == 'weibull' else 'training_diagnostics'}_{key}.json"
        )
        print(f"Fit training model: {key}", flush=True)
        try:
            model.fit(matrices[feature], outcome_labels(train, endpoint))
        except (BayesianDiagnosticsError, BayesianInitializationError, PriorPredictiveError):
            failed = run.outputs / "failed_fits" / key
            failed.mkdir(parents=True, exist_ok=True)
            failure_path = failed / "diagnostics.json"
            write_json(
                failure_path, {"diagnostics": model.diagnostics_, "workflow": model.workflow_}
            )
            artifacts = [failure_path]
            if model.idata is not None:
                model.save(failed / "posterior.nc")
                artifacts += [failed / "posterior.nc", failed / "posterior.json"]
            run.finish(artifacts, status="failed_bayesian_gate", model_key=key)
            raise
        model.save(path)
        artifacts = [path]
        if name == "weibull":
            artifacts.append(path.with_suffix(".json"))
            diagnostics = {"diagnostics": model.diagnostics_, "workflow": model.workflow_}
        elif name == "cox":
            diagnostics = {
                "proportional_hazards": model.ph_diagnostic_,
                "constant_features": model.constant_features_,
                "log_hazard_ratio": dict(
                    zip(model.active_feature_names_, model.estimator_.coef_, strict=True)
                ),
            }
        else:
            diagnostics = {"untuned": True, "n_trees": cfg.rsf.n_estimators}
        write_json(
            metric_path,
            {"endpoint": endpoint, "feature_set": feature, "model": name, **diagnostics},
        )
        artifacts.append(metric_path)
        record = {
            "endpoint": endpoint,
            "feature_set": feature,
            "model": name,
            "family": family,
            "prior_variant": prior,
            "path": str(path),
            "metrics_path": str(metric_path),
            "artifacts_sha256": {str(p): sha256_file(p) for p in artifacts},
            "source_sha256": run.metadata["source_sha256"],
        }
        index[section][key] = record
        index["artifacts_sha256"].update(record["artifacts_sha256"])
        write_json(index_path, index)
    # Compare sensitivity curves on the SAME training predictor distribution only.
    comparisons = {}
    for key, record in index["sensitivities"].items():
        reference = index["models"][f"{record['endpoint']}__weibull__{record['feature_set']}"]
        baseline = load_fit(reference)
        alternative = load_fit(record)
        X = matrices[record["feature_set"]]
        times = np.asarray(cfg.bayesian.ppc_time_grid_months, float)
        delta = alternative.predict_survival_function(
            X, times
        ) - baseline.predict_survival_function(X, times)
        comparisons[key] = {
            "split": "train",
            "times_months": times.tolist(),
            "mean_survival_difference": delta.mean(axis=0).tolist(),
            "mean_absolute_survival_difference": np.abs(delta).mean(axis=0).tolist(),
            "maximum_absolute_survival_difference": float(np.abs(delta).max()),
        }
    sensitivity_path = run.outputs / "metrics/bayesian_sensitivity.json"
    write_json(sensitivity_path, comparisons)
    index["complete"] = True
    index["artifacts_sha256"][str(sensitivity_path)] = sha256_file(sensitivity_path)
    write_json(index_path, index)
    run.finish(
        [index_path, *map(Path, index["artifacts_sha256"])], status="complete", n_fits=len(tasks)
    )
    return index_path
