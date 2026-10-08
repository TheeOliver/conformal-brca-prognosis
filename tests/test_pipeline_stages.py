"""End-to-end classical synthetic integration and one-time held-out access contract."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

from brca.calibration import run_calibration
from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT, load_config
from brca.data.eda import outcome_labels
from brca.data.load import load_cohort
from brca.data.partitions import freeze_partition_tables, load_partition
from brca.data.preprocessing import TrainingPreprocessor
from brca.data.splits import make_splits
from brca.evaluation.metrics import training_time_grid
from brca.fitting import load_fit_index
from brca.manifest import sha256_file, write_json_exclusive
from brca.models.cox import CoxPHModel
from brca.pipeline import StageRun, named_rng, write_json


def _evaluation_module():
    spec = importlib.util.spec_from_file_location(
        "stage05", PROJECT_ROOT / "scripts/05_evaluate.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def experiment(tmp_path):
    payload = yaml.safe_load(DEFAULT_CONFIG_PATH.read_text())
    payload["paths"]["processed"] = str(tmp_path / "processed")
    payload["paths"]["outputs"] = str(tmp_path / "outputs")
    payload["conformal"]["alphas"] = [0.2]
    payload["conformal"]["horizons_months"] = [120.0]
    payload["evaluation"]["bootstrap_replicates"] = 10
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(payload))
    cfg = load_config(path)
    frame = load_cohort(PROJECT_ROOT / "tests/fixtures/synthetic_metabric.csv")
    frame["event_disease_specific"] = frame.event
    frame["event_overall"] = True
    document = {
        "splits": make_splits(frame, cfg),
        "cohort_sha256": hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest(),
    }
    processed = Path(cfg.paths.processed)
    write_json_exclusive(processed / "splits.json", document)
    freeze_partition_tables(frame, document, processed / "splits.json")
    train = load_partition(cfg, "train")
    run = StageRun("synthetic_setup", path, allow_dirty=True)
    index = {
        "complete": True,
        "config_sha256": sha256_file(path),
        "split_sha256": sha256_file(processed / "splits.json"),
        "models": {},
        "preprocessors": {},
        "time_grids": {},
        "artifacts_sha256": {},
    }
    for feature in ("clinical", "clinical_molecular"):
        pre = TrainingPreprocessor(cfg, feature).fit(train)
        pre_path = run.outputs / "models" / f"{feature}.pkl"
        pre.save(pre_path)
        index["preprocessors"][feature] = str(pre_path)
        index["artifacts_sha256"][str(pre_path)] = sha256_file(pre_path)
        for endpoint in ("disease_specific", "overall"):
            y = outcome_labels(train, endpoint)
            index["time_grids"][endpoint] = training_time_grid(y, cfg).tolist()
            key = f"{endpoint}__cox__{feature}"
            model = CoxPHModel(cfg.cox, rng=named_rng(cfg.seed, key)).fit(pre.transform(train), y)
            model_path = run.outputs / "models" / f"{key}.pkl"
            model.save(model_path)
            record = {
                "endpoint": endpoint,
                "feature_set": feature,
                "model": "cox",
                "path": str(model_path),
                "artifacts_sha256": {str(model_path): sha256_file(model_path)},
            }
            index["models"][key] = record
            index["artifacts_sha256"].update(record["artifacts_sha256"])
    write_json(run.outputs / "models/fit_index.json", index)
    return path, run, index


def test_calibration_and_final_evaluation_preserve_frozen_fits(experiment, monkeypatch):
    path, run, index = experiment
    fingerprints = dict(index["artifacts_sha256"])
    calibration_path = run_calibration(run)
    calibration = json.loads(calibration_path.read_text())
    for state in calibration["states"].values():
        fitted = state["calibration"]
        if state["registered_method"] == "conservative":
            assert sum(fitted["stratum_counts"].values()) == fitted["n_cal"]
            assert set(fitted["stratum_counts"]) == {"False", "True"}
    stage = _evaluation_module()
    monkeypatch.setattr("sys.argv", ["05_evaluate.py", "--config", str(path), "--allow-dirty"])
    stage.main()
    report = json.loads((run.outputs / "metrics/performance.json").read_text())
    assert set(report) == {"disease_specific", "overall"}
    for endpoint in report.values():
        assert endpoint["bootstrap"]["sampling_design"] == "fixed_primary_event_strata"
        assert set(endpoint["paired_feature_comparisons"]) == {"cox"}
    coverage = json.loads((run.outputs / "metrics/conformal_coverage.json").read_text())
    assert len(coverage) == len(calibration["states"])
    for row in coverage.values():
        assert row["marginal"]["width_median_months"] <= 120
        assert row["marginal"]["n_cal"] == 100
    assert all(sha256_file(p) == value for p, value in fingerprints.items())
    # Recovery must use the cache even if raw held-out rows are no longer available.
    (run.processed / "partitions/test.csv").write_text("must not reopen")
    stage.main()
    assert report == json.loads((run.outputs / "metrics/performance.json").read_text())
    with pytest.raises(ValueError, match="recalibration"):
        run_calibration(run)
    prediction = next((run.outputs / "evaluation").glob("predictions__*.pkl"))
    prediction.write_bytes(prediction.read_bytes() + b"tampered")
    with pytest.raises(ValueError, match="Prediction cache changed"):
        stage.main()


def test_test_cache_rejects_protocol_change_and_tampering(experiment):
    _, run, index = experiment
    calibration_path = run_calibration(run)
    calibration = json.loads(calibration_path.read_text())
    stage = _evaluation_module()
    frame = stage.frozen_test_cache(run, index, calibration)
    assert "patient_id" not in frame
    altered = {**index, "split_sha256": "changed"}
    with pytest.raises(ValueError, match="cannot change"):
        stage.frozen_test_cache(run, altered, calibration)
    cache = run.outputs / "evaluation/test_cache.pkl"
    cache.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="cache changed"):
        stage.frozen_test_cache(run, index, calibration)
    model_path = next(iter(index["artifacts_sha256"]))
    Path(model_path).write_bytes(b"tampered")
    with pytest.raises(ValueError, match="artifact changed"):
        load_fit_index(run)
