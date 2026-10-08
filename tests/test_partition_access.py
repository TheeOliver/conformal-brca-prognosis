"""Frozen partition access and training EDA checks on synthetic rows."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from brca.config import load_config
from brca.data.eda import training_eda
from brca.data.load import load_cohort
from brca.data.partitions import freeze_partition_tables, load_partition
from brca.data.splits import make_splits
from brca.manifest import write_json_exclusive


def test_partition_access_and_tamper_detection(tmp_path):
    config = load_config()
    config.paths.processed = str(tmp_path)
    frame = load_cohort(Path(__file__).parent / "fixtures" / "synthetic_metabric.csv")
    document = {
        "splits": make_splits(frame, config),
        "cohort_sha256": hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest(),
    }
    split_path = tmp_path / "splits.json"
    write_json_exclusive(split_path, document)
    freeze_partition_tables(frame, document, split_path)
    loaded = load_partition(config, "train")
    assert loaded.patient_id.tolist() == document["splits"]["train"]
    with pytest.raises(ValueError, match="stage 05"):
        load_partition(config, "test")
    test_path = tmp_path / "partitions" / "test.csv"
    # Damaging the test file cannot affect stages loading train alone.
    test_path.write_text("unreadable")
    assert load_partition(config, "train").equals(loaded)
    with pytest.raises(ValueError, match="checksum"):
        load_partition(config, "test", allow_test=True)


def test_training_eda_excludes_patient_identifiers_and_small_subgroup_curves():
    cohort = (
        load_cohort(Path(__file__).parent / "fixtures" / "synthetic_metabric.csv").iloc[:100].copy()
    )
    cohort["event_disease_specific"] = cohort.event
    cohort["event_overall"] = cohort.event
    cohort.loc[cohort.index[:2], "pam50_subtype"] = "unknown"
    cohort.loc[cohort.index[0], "tumor_grade"] = float("nan")
    report = training_eda(cohort, load_config())
    assert report["split"] == "train" and report["n_training"] == 100
    assert "SYN-" not in json.dumps(report)
    unknown = report["survival"]["disease_specific"]["groups"]["pam50_subtype"]["unknown"]
    assert unknown == {"n": 2, "status": "suppressed_small_group"}
