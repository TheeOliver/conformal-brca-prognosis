"""Stage 02 contract checks using only synthetic cohorts and temporary outputs."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT, load_config
from brca.data.load import load_cohort
from brca.data.splits import SPLIT_NAMES, make_splits, run_split_stage
from brca.manifest import sha256_file, write_json_exclusive

FIXTURE = Path(__file__).parent / "fixtures" / "synthetic_metabric.csv"


@pytest.fixture
def cohort():
    return load_cohort(FIXTURE)


@pytest.fixture
def config():
    return load_config()


@pytest.fixture
def stage(tmp_path):
    """An isolated unborn Git repository, matching the current setup phase."""
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True, capture_output=True)
    fixture = tmp_path / "tests" / "fixtures" / FIXTURE.name
    fixture.parent.mkdir(parents=True)
    shutil.copyfile(FIXTURE, fixture)
    contents = yaml.safe_load(DEFAULT_CONFIG_PATH.read_text())
    contents["paths"]["processed"] = str(tmp_path / "processed")
    contents["paths"]["outputs"] = str(tmp_path / "artifacts")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(contents))
    return SimpleNamespace(
        root=tmp_path,
        fixture=fixture,
        config_path=config_path,
        kwargs={
            "config_path": config_path,
            "project_root": tmp_path,
            "synthetic": True,
            "allow_dirty": True,
        },
    )


def test_partitions_are_disjoint_and_exhaustive(cohort, config):
    splits = make_splits(cohort, config)
    assert tuple(splits) == SPLIT_NAMES
    flattened = [patient for ids in splits.values() for patient in ids]
    assert len(flattened) == len(set(flattened)) == len(cohort)
    assert set(flattened) == set(cohort["patient_id"])


def test_determinism_includes_row_order_and_preserves_input(cohort, config):
    original = cohort.copy(deep=True)
    expected = make_splits(cohort, config)
    assert make_splits(cohort, config) == expected
    assert make_splits(cohort.iloc[::-1], config) == expected
    assert cohort.equals(original)
    changed = deepcopy(config)
    changed.seed += 1
    assert make_splits(cohort, changed) != expected


@pytest.mark.parametrize("n", [399, 400])
def test_fraction_and_event_rate_targets(cohort, config, n):
    cohort = cohort.iloc[:n]
    labels = cohort.set_index("patient_id")["event"]
    for name, ids in make_splits(cohort, config).items():
        assert abs(len(ids) - n * getattr(config.split, name)) < 1
        assert abs(labels.loc[ids].mean() - labels.mean()) <= 1 / len(ids)
        assert labels.loc[ids].any() and not labels.loc[ids].all()


def test_non_default_config_controls_all_partition_sizes(cohort, config):
    config.split.train, config.split.calibration, config.split.test = (0.6, 0.2, 0.2)
    splits = make_splits(cohort, config)
    assert [len(splits[name]) for name in SPLIT_NAMES] == [240, 80, 80]


@pytest.mark.parametrize("rare_event", [False, True])
def test_rare_class_is_retained_when_proportional_rounding_would_omit_it(
    cohort, config, rare_event
):
    cohort = cohort.iloc[:100].copy()
    cohort["event"] = not rare_event
    cohort.loc[cohort.index[:3], "event"] = rare_event
    config.split.train, config.split.calibration, config.split.test = (0.8, 0.1, 0.1)
    labels = cohort.set_index("patient_id")["event"]
    for ids in make_splits(cohort, config).values():
        assert (labels.loc[ids] == rare_event).sum() == 1


def test_minimal_feasible_cohort_and_equal_fraction_ties(cohort, config):
    cohort = cohort.iloc[:6].copy()
    cohort["event"] = [False, True] * 3
    config.split.train = config.split.calibration = config.split.test = 1 / 3
    splits = make_splits(cohort, config)
    labels = cohort.set_index("patient_id")["event"]
    assert all(len(ids) == 2 and labels.loc[ids].sum() == 1 for ids in splits.values())


@pytest.mark.parametrize("fractions", [(0, 0.5, 0.5), (0.6, 0.3, 0.3), (np.nan, 0.5, 0.5)])
def test_invalid_fractions_rejected(cohort, config, fractions):
    config.split.train, config.split.calibration, config.split.test = fractions
    with pytest.raises(ValueError, match="fractions"):
        make_splits(cohort, config)


@pytest.mark.parametrize("seed", [-1, True, 1.5])
def test_invalid_seed_rejected(cohort, config, seed):
    config.seed = seed
    with pytest.raises(ValueError, match="seed"):
        make_splits(cohort, config)


def test_predictor_stratification_rejected(cohort, config):
    config.split.stratify_by = "pam50_subtype"
    with pytest.raises(ValueError, match="stratification"):
        make_splits(cohort, config)


@pytest.mark.parametrize("bad_id", [None, "", " ", 123])
def test_missing_or_non_string_ids_rejected(cohort, config, bad_id):
    cohort["patient_id"] = cohort["patient_id"].astype(object)
    cohort.loc[cohort.index[0], "patient_id"] = bad_id
    with pytest.raises(ValueError, match="patient_id"):
        make_splits(cohort, config)


def test_duplicate_ids_and_non_bool_events_rejected(cohort, config):
    with pytest.raises(ValueError, match="bool"):
        make_splits(cohort.assign(event=cohort.event.astype(int)), config)
    cohort.loc[cohort.index[1], "patient_id"] = cohort.loc[cohort.index[0], "patient_id"]
    with pytest.raises(ValueError, match="duplicate"):
        make_splits(cohort, config)


def test_small_infeasible_stratum_rejected(cohort, config):
    cohort["event"] = False
    cohort.loc[cohort.index[:2], "event"] = True
    with pytest.raises(ValueError, match="Too few"):
        make_splits(cohort, config)


def test_partition_too_small_for_both_classes_rejected(cohort, config):
    config.split.train, config.split.calibration, config.split.test = (0.998, 0.001, 0.001)
    with pytest.raises(ValueError, match="Too few"):
        make_splits(cohort, config)


def test_loader_preserves_numeric_looking_ids(cohort, tmp_path):
    cohort["patient_id"] = [f"{index:06d}" for index in range(len(cohort))]
    path = tmp_path / "synthetic.csv"
    cohort.to_csv(path, index=False)
    assert load_cohort(path)["patient_id"].equals(cohort["patient_id"])


def test_output_schema_and_complete_manifest(stage, monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "synthetic-test-job")
    split_path, manifest_path = run_split_stage(**stage.kwargs)
    document = json.loads(split_path.read_bytes())
    assert split_path == stage.root / "processed" / "splits.json"
    assert set(document) == {
        "schema_version",
        "seed",
        "fractions",
        "stratify_by",
        "cohort_sha256",
        "input_kind",
        "splits",
    }
    assert document["schema_version"] == 1
    assert document["seed"] == load_config(stage.config_path).seed
    assert document["stratify_by"] == "event"
    assert document["cohort_sha256"] == sha256_file(stage.fixture)
    assert set(document["splits"]) == set(SPLIT_NAMES)
    assert all(ids == sorted(ids) for ids in document["splits"].values())
    assert split_path.read_bytes().endswith(b"\n")
    manifest = json.loads(manifest_path.read_bytes())
    assert manifest_path.parent == stage.root / "artifacts" / "manifests"
    assert manifest_path.name.startswith("02_make_splits_")
    assert manifest["stage"] == "02_make_splits"
    assert manifest["git_sha"] is None and manifest["git_dirty"] is True
    assert manifest["allow_dirty"] is True
    assert manifest["config_sha256"] == sha256_file(stage.config_path)
    assert manifest["input_sha256"] == sha256_file(stage.fixture)
    assert manifest["split_sha256"] == sha256_file(split_path)
    assert manifest["input_kind"] == document["input_kind"] == "synthetic"
    assert manifest["seed"] == document["seed"]
    assert manifest["slurm_job_id"] == "synthetic-test-job"
    assert manifest["wall_clock_seconds"] >= 0
    assert manifest["timestamp_utc"].endswith("+00:00")
    assert manifest["split_reused"] is False
    assert manifest["split_sizes"] == {"train": 200, "calibration": 100, "test": 100}
    assert {"python", "numpy", "pandas", "scikit-learn", "pymc"} <= set(
        manifest["package_versions"]
    )
    assert manifest["package_versions"]["numpy"] == np.__version__


def test_frozen_reuse_is_byte_identical_without_resampling(stage, monkeypatch):
    path, first_manifest = run_split_stage(**stage.kwargs)
    original = path.read_bytes()
    original_mtime = path.stat().st_mtime_ns

    def forbidden(*args, **kwargs):
        raise AssertionError("a compatible rerun must not redraw the split")

    monkeypatch.setattr("brca.data.splits.make_splits", forbidden)
    _, second_manifest = run_split_stage(**stage.kwargs)
    assert path.read_bytes() == original
    assert path.stat().st_mtime_ns == original_mtime
    assert first_manifest != second_manifest
    assert json.loads(second_manifest.read_text())["split_reused"] is True


@pytest.mark.parametrize("change", ["seed", "fractions", "cohort"])
def test_changed_split_inputs_preserve_frozen_file(stage, change):
    path, _ = run_split_stage(**stage.kwargs)
    original = path.read_bytes()
    if change == "cohort":
        with stage.fixture.open("a") as handle:
            handle.write("\n")
    else:
        config = yaml.safe_load(stage.config_path.read_text())
        if change == "seed":
            config["seed"] += 1
        else:
            config["split"].update(train=0.6, calibration=0.2, test=0.2)
        stage.config_path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="left untouched"):
        run_split_stage(**stage.kwargs)
    assert path.read_bytes() == original
    assert len(list((stage.root / "artifacts" / "manifests").glob("*.json"))) == 1


@pytest.mark.parametrize("corruption", ["invalid_json", "overlap", "missing", "unknown_id"])
def test_invalid_frozen_artifact_is_not_repaired_silently(stage, corruption):
    path, _ = run_split_stage(**stage.kwargs)
    document = json.loads(path.read_text())
    if corruption == "invalid_json":
        path.write_text("{")
    else:
        splits = document["splits"]
        if corruption == "overlap":
            splits["test"][0] = splits["train"][0]
            splits["test"].sort()
        elif corruption == "missing":
            splits["test"].pop()
        else:
            splits["test"][0] = "unknown-synthetic-id"
            splits["test"].sort()
        path.write_text(json.dumps(document))
    before = path.read_bytes()
    with pytest.raises(ValueError, match="incompatible or invalid") as exc:
        run_split_stage(**stage.kwargs)
    assert "unknown-synthetic-id" not in str(exc.value)
    assert path.read_bytes() == before


def test_dirty_tree_requires_explicit_override_before_writing(stage):
    with pytest.raises(ValueError, match="--allow-dirty"):
        run_split_stage(**(stage.kwargs | {"allow_dirty": False}))
    assert not (stage.root / "processed").exists()
    assert not (stage.root / "artifacts").exists()


def test_missing_processed_input_never_falls_back_to_synthetic(stage):
    with pytest.raises(ValueError, match="Stage 01"):
        run_split_stage(**(stage.kwargs | {"synthetic": False}))
    assert not (stage.root / "processed").exists()


def test_explicit_synthetic_mode_keeps_label_for_symlinked_fixture(stage):
    target = stage.root / "synthetic.csv"
    stage.fixture.rename(target)
    stage.fixture.symlink_to(target)
    split_path, manifest_path = run_split_stage(**stage.kwargs)
    assert json.loads(split_path.read_text())["input_kind"] == "synthetic"
    assert json.loads(manifest_path.read_text())["input_kind"] == "synthetic"


def test_atomic_writer_never_clobbers_existing_file(tmp_path):
    path = tmp_path / "frozen.json"
    write_json_exclusive(path, {"original": True})
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        write_json_exclusive(path, {"original": False})
    assert path.read_bytes() == before
    assert list(tmp_path.iterdir()) == [path]


def test_cli_synthetic_run(stage):
    completed = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "02_make_splits.py"),
            "--config",
            str(stage.config_path),
            "--synthetic",
            "--allow-dirty",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "Frozen split:" in completed.stdout and "Run manifest:" in completed.stdout
    assert "SYN-" not in completed.stdout
    assert (stage.root / "processed" / "splits.json").is_file()
