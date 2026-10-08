"""Synthetic tampering checks for the read-only provenance audit."""

from __future__ import annotations

import importlib.util
import io
import tarfile
from pathlib import Path

from brca.manifest import sha256_file

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/audit_results.py"
spec = importlib.util.spec_from_file_location("result_audit", SCRIPT)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_digest_detects_artifact_change_and_missing_file(tmp_path):
    artifact = tmp_path / "aggregate.json"
    artifact.write_text('{"n": 100}')
    expected = sha256_file(artifact)
    checks = audit.IntegrityChecks()
    assert checks.digest(artifact, expected, "original")
    artifact.write_text('{"n": 101}')
    assert not checks.digest(artifact, expected, "edited")
    artifact.unlink()
    assert not checks.digest(artifact, expected, "missing")
    assert checks.failures == ["edited", "missing"]


def test_manifest_requires_current_config_and_all_output_digests():
    original = {
        "stage": "04_conformal",
        "status": "complete",
        "config_sha256": "config",
        "artifacts_sha256": {"calibration.json": "original"},
    }
    manifests = [(Path("04_a.json"), original)]
    checks = audit.IntegrityChecks()
    assert audit.select_manifest(
        manifests, "04_conformal", "config", {"calibration.json": "original"}, checks
    )
    assert (
        audit.select_manifest(
            manifests, "04_conformal", "config", {"calibration.json": "edited"}, checks
        )
        is None
    )
    assert (
        audit.select_manifest(
            manifests, "04_conformal", "changed_config", {"calibration.json": "original"}, checks
        )
        is None
    )
    assert len(checks.failures) == 2


def test_source_archive_rejects_patient_root_and_unsafe_paths(tmp_path):
    path = tmp_path / "source.tar.gz"
    content = b"synthetic fixture, not patient data"
    with tarfile.open(path, "w:gz") as archive:
        for name in ("src/brca/example.py", "data/forbidden.txt", "../escape.py"):
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    checks = audit.IntegrityChecks()
    hashes = audit.archive_hashes(path, sha256_file(path), checks)
    assert list(hashes) == ["src/brca/example.py"]
    assert len(checks.failures) == 2
    assert not (tmp_path.parent / "escape.py").exists()
