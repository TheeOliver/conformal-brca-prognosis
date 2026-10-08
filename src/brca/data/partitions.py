"""Physically separate cohort partitions so downstream stages never load test rows."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

from brca.config import PROJECT_ROOT
from brca.data.load import load_cohort
from brca.manifest import sha256_file, write_json_exclusive


def archive_synthetic_split(processed: Path, outputs: Path) -> Path | None:
    """Preserve the development split before the explicitly requested real transition."""
    path = processed / "splits.json"
    if not path.exists():
        return None
    document = json.loads(path.read_text())
    if document.get("input_kind") != "synthetic":
        return None
    archive = processed / "archive" / f"synthetic_{sha256_file(path)}"
    archive.mkdir(parents=True, exist_ok=False)
    shutil.move(str(path), str(archive / "splits.json"))
    partitions = processed / "partitions"
    if partitions.exists():
        shutil.move(str(partitions), str(archive / "partitions"))
    for manifest in (outputs / "manifests").glob("02_make_splits_*.json"):
        if json.loads(manifest.read_text()).get("input_kind") == "synthetic":
            shutil.move(str(manifest), str(archive / manifest.name))
    return archive


def freeze_partition_tables(cohort: pd.DataFrame, document: dict, split_path: Path) -> dict:
    """Called only by stage 02, immediately after the split is drawn/validated."""
    directory = split_path.parent / "partitions"
    directory.mkdir(parents=True, exist_ok=True)
    indexed = cohort.set_index("patient_id", drop=False)
    hashes = {}
    for name, ids in document["splits"].items():
        content = indexed.loc[ids].to_csv(index=False, lineterminator="\n").encode()
        path = directory / f"{name}.csv"
        if path.exists() and path.read_bytes() != content:
            raise ValueError("Frozen partition table differs from the registered split.")
        if not path.exists():
            with path.open("xb") as handle:
                handle.write(content)
        hashes[name] = hashlib.sha256(content).hexdigest()
    payload = {
        "cohort_sha256": document["cohort_sha256"],
        "split_sha256": sha256_file(split_path),
        "partitions_sha256": hashes,
    }
    metadata = directory / "checksums.json"
    if metadata.exists():
        if json.loads(metadata.read_text()) != payload:
            raise ValueError("Frozen partition checksums are inconsistent.")
    else:
        write_json_exclusive(metadata, payload)
    return payload


def load_partition(config, name: str, *, allow_test: bool = False) -> pd.DataFrame:
    """Load one frozen table; only final stage 05 may request ``allow_test=True``."""
    if name not in {"train", "calibration", "test"} or (name == "test" and not allow_test):
        raise ValueError("Test patients may only be loaded by final stage 05.")
    root = PROJECT_ROOT / config.paths.processed
    split_path = root / "splits.json"
    split = json.loads(split_path.read_text())
    checksums = json.loads((root / "partitions" / "checksums.json").read_text())
    path = root / "partitions" / f"{name}.csv"
    if checksums["split_sha256"] != sha256_file(split_path) or checksums["partitions_sha256"][
        name
    ] != sha256_file(path):
        raise ValueError("Frozen split or partition checksum mismatch.")
    cohort = load_cohort(path)
    if cohort["patient_id"].tolist() != split["splits"][name]:
        raise ValueError("Partition patient membership differs from the frozen split.")
    return cohort
