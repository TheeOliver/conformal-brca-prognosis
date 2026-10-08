"""Calibration-only conformal orchestration; models and preprocessing stay frozen."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from brca.conformal.survival import calibrate_intervals, calibrate_ipcw_intervals
from brca.data.eda import outcome_labels
from brca.data.partitions import load_partition
from brca.data.preprocessing import TrainingPreprocessor
from brca.fitting import load_fit, load_fit_index
from brca.manifest import sha256_file
from brca.pipeline import StageRun, write_json


def interval_key(model_key: str, alpha: float, horizon: float, method: str) -> str:
    return f"{model_key}__alpha{alpha:g}__horizon{horizon:g}__{method}"


def model_quantiles(model, X, alphas) -> dict[float, np.ndarray]:
    probabilities = [p for alpha in alphas for p in (alpha / 2, 1 - alpha / 2)]
    values = model.predict_quantiles(X, np.asarray(probabilities))
    return {alpha: values[:, 2 * j : 2 * j + 2] for j, alpha in enumerate(alphas)}


def run_calibration(run: StageRun) -> Path:
    if (run.outputs / "evaluation/plan.json").exists():
        raise ValueError("Test evaluation has opened; recalibration is prohibited.")
    cfg = run.config
    index = load_fit_index(run)
    train = load_partition(cfg, "train")
    calibration = load_partition(cfg, "calibration")
    states = {}
    for key, record in index["models"].items():
        model = load_fit(record)
        preprocessor = TrainingPreprocessor.load(
            Path(index["preprocessors"][record["feature_set"]])
        )
        X = preprocessor.transform(calibration)
        y_train = outcome_labels(train, record["endpoint"])
        y_cal = outcome_labels(calibration, record["endpoint"])
        quantiles = model_quantiles(model, X, cfg.conformal.alphas)
        for alpha, base in quantiles.items():
            for horizon in cfg.conformal.horizons_months:
                lower, upper = np.minimum(base, horizon).T
                for method in cfg.conformal.methods:
                    if method == "conservative":
                        fitted = calibrate_intervals(
                            y_cal,
                            lower,
                            upper,
                            alpha,
                            horizon=horizon,
                            strata=calibration["event"].to_numpy(bool),
                            expected_strata=(False, True),
                        )
                    elif method == "ipcw_sensitivity":
                        fitted = calibrate_ipcw_intervals(
                            y_cal,
                            lower,
                            upper,
                            y_train,
                            alpha,
                            horizon=horizon,
                            min_censoring_survival=cfg.conformal.min_censoring_survival,
                        )
                    else:
                        raise ValueError("Unregistered conformal method.")
                    states[interval_key(key, alpha, horizon, method)] = {
                        "model_key": key,
                        "endpoint": record["endpoint"],
                        "feature_set": record["feature_set"],
                        "model": record["model"],
                        "registered_method": method,
                        "calibration": fitted.to_dict(),
                    }
        print(f"Calibrated frozen fit: {key}", flush=True)
    path = run.outputs / "models/conformal_calibration.json"
    write_json(
        path,
        {
            "config_sha256": run.metadata["config_sha256"],
            "split_sha256": index["split_sha256"],
            "fit_index_sha256": sha256_file(run.outputs / "models/fit_index.json"),
            "calibration_partition_sha256": sha256_file(
                run.processed / "partitions/calibration.csv"
            ),
            "states": states,
        },
    )
    run.finish([path], status="complete", n_calibrations=len(states), test_accessed=False)
    return path


def load_calibration(run: StageRun, index: dict) -> dict:
    path = run.outputs / "models/conformal_calibration.json"
    state = json.loads(path.read_text())
    if (
        state["config_sha256"] != run.metadata["config_sha256"]
        or state["split_sha256"] != index["split_sha256"]
        or state["fit_index_sha256"] != sha256_file(run.outputs / "models/fit_index.json")
    ):
        raise ValueError("Calibration belongs to a different frozen analysis.")
    return state
