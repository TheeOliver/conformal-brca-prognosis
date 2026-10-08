"""Shared pipeline provenance and deterministic output helpers."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import tarfile
from datetime import UTC, datetime
from pathlib import Path
from tempfile import NamedTemporaryFile
from time import perf_counter

import numpy as np

from brca.config import PROJECT_ROOT, load_config
from brca.manifest import git_state, package_versions, sha256_file, write_json_exclusive


def json_safe(value):
    """JSON has no NaN or infinity; undefined scalar summaries become null."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path: Path, value: dict) -> None:
    """Atomically replace a regenerable aggregate output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(json_safe(value), handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)


def named_rng(seed: int, *names: str) -> np.random.Generator:
    """Stable independent streams: adding a model never shifts another model's RNG."""
    digest = hashlib.sha256("/".join(names).encode()).digest()
    entropy = np.frombuffer(digest, dtype="<u4").tolist()
    return np.random.default_rng(np.random.SeedSequence([seed, *entropy]))


def _source_archive(root: Path, directory: Path) -> tuple[str, Path]:
    """Preserve runnable source/config bytes even for explicitly allowed dirty runs.

    Only named code/config roots enter the archive; no patient data or outputs.
    """
    paths = []
    for name in ("src", "scripts", "config"):
        paths.extend(
            path for path in (root / name).rglob("*") if path.suffix in {".py", ".yaml", ".sh"}
        )
    paths.extend(root / name for name in ("pyproject.toml", "uv.lock") if (root / name).is_file())
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        for path in sorted(paths):
            raw = path.read_bytes()
            info = tarfile.TarInfo(str(path.relative_to(root)))
            info.size = len(raw)
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(raw))
    content = gzip.compress(buffer.getvalue(), mtime=0)
    digest = hashlib.sha256(content).hexdigest()
    target = directory / f"source_{digest}.tar.gz"
    if not target.exists():
        directory.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as handle:
            handle.write(content)
    return digest, target


class StageRun:
    """Record config, source snapshot, inputs, outputs and SLURM execution."""

    def __init__(self, stage: str, config_path: Path | str, *, allow_dirty: bool = False):
        self.started = perf_counter()
        self.stage = stage
        self.config_path = Path(config_path).resolve()
        self.config = load_config(self.config_path)
        self.outputs = PROJECT_ROOT / self.config.paths.outputs
        self.processed = PROJECT_ROOT / self.config.paths.processed
        self.timestamp = datetime.now(UTC)
        state = git_state(PROJECT_ROOT)
        if state["git_dirty"] and not allow_dirty:
            raise ValueError("Uncommitted tree: pass --allow-dirty to record a source snapshot.")
        source_hash, source_path = _source_archive(PROJECT_ROOT, self.outputs / "manifests")
        self.metadata = {
            "schema_version": 1,
            "stage": stage,
            "timestamp_utc": self.timestamp.isoformat(),
            **state,
            "allow_dirty": allow_dirty,
            "source_sha256": source_hash,
            "source_archive": str(source_path),
            "config_sha256": sha256_file(self.config_path),
            "config_path": str(self.config_path),
            "seed": self.config.seed,
            "package_versions": package_versions(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        }

    def finish(self, artifacts: list[Path], **metadata) -> Path:
        manifest = (
            self.metadata
            | metadata
            | {
                "wall_clock_seconds": perf_counter() - self.started,
                "artifacts_sha256": {str(path): sha256_file(path) for path in artifacts},
            }
        )
        path = (
            self.outputs
            / "manifests"
            / (f"{self.stage}_{self.timestamp.strftime('%Y%m%dT%H%M%S%fZ')}.json")
        )
        write_json_exclusive(path, json_safe(manifest))
        return path
