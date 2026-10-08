"""Synthetic integration and immutable provenance for the exploratory experiment."""

from __future__ import annotations

import json

import numpy as np
import pytest
import yaml

from brca.calibration import run_calibration
from brca.manifest import sha256_file
from brca.pipeline import StageRun, write_json
from brca.reliability import (
    freeze_reliability_plan,
    prepare_reliability,
    training_support,
)
from tests.test_pipeline_stages import _evaluation_module, experiment  # noqa: F401


def test_support_audit_handles_horizon_without_censoring_support():
    y = np.array([(True, 1.0), (True, 2.0), (False, 3.0)], dtype=[("event", "?"), ("time", "f8")])
    summary = training_support(y, [2.0, 4.0], 0.05)
    assert summary["horizons"]["2.0"]["ipcw_support_status"] == "supported"
    beyond = summary["horizons"]["4.0"]
    assert beyond["ipcw_support_status"] == "unsupported"
    assert beyond["restricted_known_target_effective_sample_size"] is None
    assert beyond["weighted_horizon_mass"] is None


def test_reliability_pipeline_and_provenance(experiment, monkeypatch, tmp_path):  # noqa: F811
    path, run, index = experiment
    index["training_sha256"] = sha256_file(run.processed / "partitions/train.csv")
    write_json(run.outputs / "models/fit_index.json", index)
    run_calibration(run)
    stage = _evaluation_module()
    monkeypatch.setattr("sys.argv", ["05_evaluate.py", "--config", str(path), "--allow-dirty"])
    stage.main()
    qin_settings = tmp_path / "qin.yaml"
    qin_settings.write_text(
        yaml.safe_dump(
            {
                "bootstrap_replicates": 10,
                "minimum_bootstrap_success_fraction": 0.5,
                "minimum_censoring_survival_for_event_pool": 0.0,
            }
        )
    )
    stage.run_qin_exploratory(
        StageRun("synthetic_qin", path, allow_dirty=True), index, qin_settings
    )
    settings_path = tmp_path / "reliability.yaml"
    settings = {
        "ridge_alpha": 1.0,
        "ridge_bootstrap_replicates": 10,
        "diagnostic_bootstrap_replicates": 10,
        "report_alpha": 0.2,
    }
    settings_path.write_text(yaml.safe_dump(settings))
    reference_paths = [
        run.outputs / "evaluation/plan.json",
        run.outputs / "qin_exploratory/coverage.json",
        run.outputs / "metrics/conformal_coverage.json",
    ]
    reference_sha = {p: sha256_file(p) for p in reference_paths}
    monkeypatch.setattr(
        "sys.argv",
        [
            "05_evaluate.py",
            "--config",
            str(path),
            "--allow-dirty",
            "--reliability-exploratory",
            "--reliability-config",
            str(settings_path),
        ],
    )
    stage.main()
    directory = run.outputs / "uncertainty_reliability"
    coverage = json.loads((directory / "coverage.json").read_text())
    assert len(coverage) == 8
    assert {row["method"] for row in coverage.values()} == {
        "conservative_lower_bound",
        "qin_bootstrap_ridge_cox",
    }
    for row in coverage.values():
        if row["method"] == "conservative_lower_bound":
            assert row["marginal"]["n_cal"] == 100
        else:
            assert row["bootstrap_diagnostics"]["bootstrap_successful"] == 10
    assert all(sha256_file(p) == digest for p, digest in reference_sha.items())
    prepared = prepare_reliability(run, index, settings)
    assert prepared["n_cal"] == 100
    # A changed protocol must be rejected before any new test access.
    settings_path.write_text(yaml.safe_dump({**settings, "ridge_alpha": 2.0}))
    with pytest.raises(ValueError, match="protocol changed"):
        freeze_reliability_plan(run, index, settings_path)
    settings_path.write_text(yaml.safe_dump(settings))
    model_path = next(iter(prepared["artifacts_sha256"]))
    with open(model_path, "ab") as handle:
        handle.write(b"tampered")
    with pytest.raises(ValueError, match="preparation cache changed"):
        prepare_reliability(run, index, settings)
    calibration = run.processed / "partitions/calibration.csv"
    calibration.write_text("tampered")
    with pytest.raises(ValueError, match="inputs differ"):
        freeze_reliability_plan(run, index, settings_path)
