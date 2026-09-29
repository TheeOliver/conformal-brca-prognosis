"""The single event-stratified partition, frozen as patient IDs by stage 02.

Only patient IDs and the event indicator determine assignments. No predictor
transforms or outcome analyses belong here.
"""

from __future__ import annotations

import json
import os
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from brca.config import PROJECT_ROOT, load_config
from brca.data.load import load_cohort
from brca.data.schema import validate
from brca.manifest import git_state, package_versions, sha256_file, write_json_exclusive

SPLIT_NAMES = ("train", "calibration", "test")
SCHEMA_VERSION = 1


def _protocol(config: SimpleNamespace) -> dict[str, Any]:
    seed = config.seed
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("config.seed must be a non-negative integer.")
    if config.split.stratify_by != "event" or config.outcome.event_column != "event":
        raise ValueError("Stage 02 requires stratification by the bool event indicator.")
    fractions = {name: float(getattr(config.split, name)) for name in SPLIT_NAMES}
    values = np.array(list(fractions.values()))
    if not np.isfinite(values).all() or (values <= 0).any() or (values >= 1).any():
        raise ValueError("Split fractions must be finite and strictly between zero and one.")
    if not np.isclose(values.sum(), 1, rtol=0, atol=np.finfo(float).eps * len(values)):
        raise ValueError("Split fractions must sum to one.")
    return {"seed": seed, "fractions": fractions, "stratify_by": config.split.stratify_by}


def _allocation(cohort: pd.DataFrame, protocol: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    validate(cohort)
    ids = cohort["patient_id"]
    if not ids.map(lambda value: isinstance(value, str) and bool(value.strip())).all():
        raise ValueError("patient_id must contain non-empty strings, with no missing values.")

    targets = len(cohort) * np.array([protocol["fractions"][name] for name in SPLIT_NAMES])
    sizes = np.floor(targets).astype(int)
    # Largest remainder rounding keeps the total exact; ties follow SPLIT_NAMES.
    priority = np.argsort(-(targets - sizes), kind="stable")
    sizes[priority[: len(cohort) - sizes.sum()]] += 1
    n_events = int(cohort["event"].sum())
    if (sizes < 2).any() or min(n_events, len(cohort) - n_events) < len(SPLIT_NAMES):
        raise ValueError("Too few patients to place events and censored rows in every split.")

    # Bounded apportionment also handles small, rare-event strata without dropping
    # rows or silently returning a partition containing just one outcome class.
    event_targets = n_events * sizes / len(cohort)
    events = np.clip(np.floor(event_targets).astype(int), 1, sizes - 1)
    while events.sum() < n_events:
        deficit = np.where(events < sizes - 1, event_targets - events, -np.inf)
        events[np.argmax(deficit)] += 1
    while events.sum() > n_events:
        surplus = np.where(events > 1, events - event_targets, -np.inf)
        events[np.argmax(surplus)] -= 1
    return sizes, events


def make_splits(cohort: pd.DataFrame, config: SimpleNamespace) -> dict[str, list[str]]:
    """Draw deterministic partitions, invariant to the input row order.

    Partition sizes use largest-remainder rounding. Each receives both event
    classes; class counts follow proportional allocation subject to that constraint.
    Tiny cohorts for which that is impossible fail before any file is written.
    """
    protocol = _protocol(config)
    sizes, events = _allocation(cohort, protocol)
    rng = np.random.default_rng(protocol["seed"])
    splits: dict[str, list[str]] = {name: [] for name in SPLIT_NAMES}
    for observed, counts in ((False, sizes - events), (True, events)):
        ids = sorted(cohort.loc[cohort["event"] == observed, "patient_id"])
        shuffled = rng.permutation(ids)
        start = 0
        for name, count in zip(SPLIT_NAMES, counts, strict=True):
            splits[name].extend(shuffled[start : start + count].tolist())
            start += count
    return {name: sorted(ids) for name, ids in splits.items()}


def _read_frozen(path: Path, metadata: dict[str, Any], cohort: pd.DataFrame) -> dict[str, Any]:
    """Validate and reuse a saved partition without drawing another one."""
    try:
        document = json.loads(path.read_text())
        if not isinstance(document, dict) or set(document) != {*metadata, "splits"}:
            raise ValueError
        if any(document[key] != value for key, value in metadata.items()):
            raise ValueError
        splits = document["splits"]
        if not isinstance(splits, dict) or set(splits) != set(SPLIT_NAMES):
            raise ValueError
        all_ids = []
        sizes, events = _allocation(cohort, metadata)
        labels = cohort.set_index("patient_id")["event"]
        for index, name in enumerate(SPLIT_NAMES):
            ids = splits[name]
            if not isinstance(ids, list) or any(not isinstance(value, str) for value in ids):
                raise ValueError
            if ids != sorted(ids) or len(ids) != sizes[index]:
                raise ValueError
            if labels.loc[ids].sum() != events[index]:
                raise ValueError
            all_ids.extend(ids)
        if len(set(all_ids)) != len(all_ids) or set(all_ids) != set(labels.index):
            raise ValueError
    except (ValueError, KeyError, TypeError):
        raise ValueError(
            "Existing splits.json is incompatible or invalid; it was left untouched. "
            "Use a separate processed/output directory for a new experiment."
        ) from None
    return document


def run_split_stage(
    config_path: Path | str,
    *,
    input_path: Path | str | None = None,
    synthetic: bool = False,
    allow_dirty: bool = False,
    project_root: Path = PROJECT_ROOT,
) -> tuple[Path, Path]:
    """Freeze the split and write provenance; return split and manifest paths.

    Relative configured paths are rooted at the project, as in stage 00. Synthetic
    input is explicit; the default expects stage 01's processed ``cohort.csv``.
    """
    started = perf_counter()
    config_path = Path(config_path).resolve()
    config_hash = sha256_file(config_path)
    config = load_config(config_path)
    protocol = _protocol(config)
    state = git_state(project_root)
    if state["git_dirty"] and not allow_dirty:
        raise ValueError(
            "Working tree is dirty or has no commits; --allow-dirty is required "
            "for a development run whose outputs are not thesis results."
        )
    if synthetic and input_path is not None:
        raise ValueError("Choose --synthetic or --input, not both.")
    fixture_root = project_root / "tests" / "fixtures"
    if synthetic:
        source = fixture_root / "synthetic_metabric.csv"
    elif input_path is not None:
        source = Path(input_path).resolve()
    else:
        source = project_root / config.paths.processed / "cohort.csv"
    source = source.resolve()
    if not source.is_file():
        raise ValueError(
            "Cohort input is missing. Stage 01 must produce processed/cohort.csv; "
            "for a synthetic check generate the fixture and pass --synthetic."
        )
    input_kind = (
        "synthetic" if synthetic or source.is_relative_to(fixture_root.resolve()) else "processed"
    )
    cohort_hash = sha256_file(source)
    try:
        cohort = load_cohort(source)
    except (pd.errors.ParserError, UnicodeError):
        raise ValueError("Cohort CSV could not be parsed; no patient rows were logged.") from None
    # Detect an input edited between hashing and loading instead of recording false provenance.
    if sha256_file(source) != cohort_hash or sha256_file(config_path) != config_hash:
        raise ValueError(
            "Input or config changed while stage 02 was reading it; no output written."
        )
    metadata = {
        "schema_version": SCHEMA_VERSION,
        **protocol,
        "cohort_sha256": cohort_hash,
        "input_kind": input_kind,
    }
    split_path = project_root / config.paths.processed / "splits.json"
    if split_path.exists():
        document = _read_frozen(split_path, metadata, cohort)
        reused = True
    else:
        document = {**metadata, "splits": make_splits(cohort, config)}
        try:
            write_json_exclusive(split_path, document)
            reused = False
        except FileExistsError:
            document = _read_frozen(split_path, metadata, cohort)
            reused = True
    timestamp = datetime.now(UTC)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "stage": "02_make_splits",
        "timestamp_utc": timestamp.isoformat(),
        **state,
        "allow_dirty": allow_dirty,
        "config_sha256": config_hash,
        "config_path": str(config_path),
        "package_versions": package_versions(),
        "seed": config.seed,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "input_kind": input_kind,
        "input_path": str(source),
        "input_sha256": cohort_hash,
        "split_path": str(split_path.resolve()),
        "split_sha256": sha256_file(split_path),
        "split_reused": reused,
        "split_sizes": {name: len(ids) for name, ids in document["splits"].items()},
        "wall_clock_seconds": perf_counter() - started,
    }
    manifest_path = (
        project_root
        / config.paths.outputs
        / "manifests"
        / f"02_make_splits_{timestamp.strftime('%Y%m%dT%H%M%S%fZ')}.json"
    )
    write_json_exclusive(manifest_path, manifest)
    return split_path, manifest_path
