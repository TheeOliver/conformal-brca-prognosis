"""Small provenance helpers shared by pipeline stages."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from importlib.metadata import distributions
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any


def sha256_file(path: Path | str) -> str:
    """Hash the exact bytes used by a stage."""
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def git_state(project_root: Path) -> dict[str, Any]:
    """Record an unborn HEAD honestly; fail if Git provenance is unavailable."""

    def git(*args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", "-C", str(project_root), *args],
            capture_output=True,
            text=True,
            check=False,
        )

    status = git("status", "--porcelain", "--untracked-files=normal")
    if status.returncode:
        raise RuntimeError("Cannot establish Git provenance for this run.")
    head = git("rev-parse", "--verify", "HEAD")
    return {
        "git_sha": head.stdout.strip() if head.returncode == 0 else None,
        "git_dirty": bool(status.stdout.strip()) or head.returncode != 0,
    }


def package_versions() -> dict[str, str]:
    """Record the resolved environment without importing expensive libraries."""
    versions = {"python": platform.python_version()}
    versions.update({dist.metadata["Name"]: dist.version for dist in distributions()})
    return dict(sorted(versions.items()))


def write_json_exclusive(path: Path, payload: dict[str, Any]) -> None:
    """Publish complete JSON atomically, refusing to replace an existing file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(content.encode("utf-8"))
            handle.flush()
            os.fsync(handle.fileno())
            os.link(temporary, path)
        finally:
            temporary.unlink()
