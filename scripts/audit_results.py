"""Audit saved result provenance without decoding patient or model objects.

Run through SLURM after stage 06 with ``--require-exports``. This verifies artifact
integrity and source/environment bindings; it does not rerun statistical analyses.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tarfile
import tomllib
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from brca.compute_guard import LoginNodeError, require_compute_node
from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT, load_config
from brca.manifest import sha256_file


def normalize_package(name: str) -> str:
    return name.lower().replace("_", "-").replace(".", "-")


class IntegrityChecks:
    """Collect all failed checks rather than stopping at the first mismatch."""

    def __init__(self):
        self.count = 0
        self.failures = []

    def check(self, condition: bool, label: str) -> bool:
        self.count += 1
        if not condition:
            self.failures.append(label)
        return condition

    def digest(self, path: Path | str, expected: str, label: str) -> bool:
        path = Path(path)
        return self.check(path.is_file() and sha256_file(path) == expected, label)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def select_manifest(
    manifests: list[tuple[Path, dict]],
    stage: str,
    config_hash: str,
    required_artifacts: dict[str, str],
    checks: IntegrityChecks,
) -> tuple[Path, dict] | None:
    """Select the newest manifest that binds every required current artifact."""
    candidates = [
        (path, manifest)
        for path, manifest in manifests
        if manifest.get("stage") == stage
        and manifest.get("config_sha256") == config_hash
        and (stage == "05_evaluate" or manifest.get("status") == "complete")
        and all(
            manifest.get("artifacts_sha256", {}).get(path) == digest
            for path, digest in required_artifacts.items()
        )
    ]
    checks.check(bool(candidates), f"matching producing manifest: {stage}")
    return max(candidates, key=lambda item: item[0].name) if candidates else None


def archive_hashes(path: Path, digest: str, checks: IntegrityChecks) -> dict[str, str]:
    """Hash archived source bytes in memory without extracting archive entries."""
    if not checks.digest(path, digest, f"source archive digest: {path.name}"):
        return {}
    hashes = {}
    allowed_roots = {"src", "scripts", "config", "pyproject.toml", "uv.lock"}
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            name = Path(member.name)
            safe = (
                bool(name.parts)
                and not name.is_absolute()
                and ".." not in name.parts
                and name.parts[0] in allowed_roots
                and member.isfile()
                and member.name not in hashes
            )
            if checks.check(safe, f"source archive member allowed: {member.name}"):
                with archive.extractfile(member) as stream:
                    hashes[member.name] = hashlib.file_digest(stream, "sha256").hexdigest()
    return hashes


def _matching_stage(manifests, matched, checks, stage, config_hash, artifacts, *, required=True):
    if not required and not any(m.get("stage") == stage for _, m in manifests):
        return
    existing = all(path.is_file() for path in artifacts)
    checks.check(existing, f"required stage outputs exist: {stage}")
    if existing:
        candidate = select_manifest(
            manifests,
            stage,
            config_hash,
            {str(path): sha256_file(path) for path in artifacts},
            checks,
        )
        if candidate:
            matched[stage] = candidate


def audit_results(root: Path, config_path: Path, *, require_exports: bool = False) -> dict:
    """Audit aggregate metadata and byte digests only; never deserialize patient rows."""
    checks = IntegrityChecks()
    cfg = load_config(config_path)
    outputs = root / cfg.paths.outputs
    processed = root / cfg.paths.processed
    index_path = outputs / "models/fit_index.json"
    calibration_path = outputs / "models/conformal_calibration.json"
    index, calibration = read_json(index_path), read_json(calibration_path)
    partition_metadata = read_json(processed / "partitions/checksums.json")
    config_hash = sha256_file(config_path)
    lock_bytes = (root / "uv.lock").read_bytes()
    lock_hash = hashlib.sha256(lock_bytes).hexdigest()
    locked = {}
    for package in tomllib.loads(lock_bytes.decode())["package"]:
        locked.setdefault(normalize_package(package["name"]), set()).add(package["version"])

    checks.check(index["complete"], "fit index complete")
    checks.check(
        index["config_sha256"] == calibration["config_sha256"] == config_hash,
        "fit and calibration config bindings",
    )
    checks.check(
        index["split_sha256"] == calibration["split_sha256"] == partition_metadata["split_sha256"],
        "split metadata bindings",
    )
    checks.check(
        index["training_sha256"] == partition_metadata["partitions_sha256"]["train"],
        "training partition metadata binding",
    )
    checks.check(
        calibration["calibration_partition_sha256"]
        == partition_metadata["partitions_sha256"]["calibration"],
        "calibration partition metadata binding",
    )
    checks.digest(index_path, calibration["fit_index_sha256"], "calibration fit-index binding")
    endpoints = [cfg.outcome.primary_endpoint, *cfg.outcome.sensitivity_endpoints]
    features = list(vars(cfg.feature_sets))
    expected_models = {
        (endpoint, feature, model)
        for endpoint in endpoints
        for feature in features
        for model in ("cox", "rsf", "weibull")
    }
    checks.check(
        {(r["endpoint"], r["feature_set"], r["model"]) for r in index["models"].values()}
        == expected_models
        and len(index["models"]) == len(expected_models),
        "registered model matrix",
    )
    expected_sensitivities = {
        (cfg.outcome.primary_endpoint, feature, "weibull", family, prior)
        for feature in features
        for family, prior in (("weibull", "sensitivity"), ("lognormal", "default"))
    }
    checks.check(
        {
            (r["endpoint"], r["feature_set"], r["model"], r["family"], r["prior_variant"])
            for r in index["sensitivities"].values()
        }
        == expected_sensitivities
        and len(index["sensitivities"]) == len(expected_sensitivities),
        "registered sensitivity matrix",
    )
    expected_calibrations = {
        (key, alpha, horizon, method)
        for key in index["models"]
        for alpha in cfg.conformal.alphas
        for horizon in cfg.conformal.horizons_months
        for method in cfg.conformal.methods
    }
    checks.check(
        {
            (
                r["model_key"],
                r["calibration"]["alpha"],
                r["calibration"]["horizon"],
                r["registered_method"],
            )
            for r in calibration["states"].values()
        }
        == expected_calibrations
        and len(calibration["states"]) == len(expected_calibrations),
        "registered calibration matrix",
    )
    for path, expected in index["artifacts_sha256"].items():
        checks.digest(path, expected, f"registered fit artifact: {Path(path).name}")
    for path, expected in index["fit_source_sha256"].items():
        checks.digest(path, expected, f"current fitting source: {Path(path).name}")

    manifests = [(p, read_json(p)) for p in (outputs / "manifests").glob("*.json")]
    preparation_matches = [
        (path, manifest)
        for path, manifest in manifests
        if manifest.get("stage") == "01_prepare_data"
        and manifest.get("status", "complete") == "complete"
        and manifest.get("config_sha256") == config_hash
        and manifest.get("input_kind") == "real_metabric"
        and manifest.get("artifacts_sha256", {}).get(str(processed / "cohort.csv"))
        == partition_metadata["cohort_sha256"]
    ]
    split_matches = [
        (path, manifest)
        for path, manifest in manifests
        if manifest.get("stage") == "02_make_splits"
        and manifest.get("config_sha256") == config_hash
        and manifest.get("input_kind") == "processed"
        and manifest.get("input_sha256") == partition_metadata["cohort_sha256"]
        and manifest.get("partition_metadata") == partition_metadata
        and manifest.get("split_sha256") == index["split_sha256"]
    ]
    checks.check(bool(preparation_matches), "real METABRIC preparation metadata binding")
    checks.check(bool(split_matches), "frozen partition metadata from processed real cohort")
    data_origin = {}
    for label, candidates in (("preparation", preparation_matches), ("split", split_matches)):
        if candidates:
            path, manifest = max(candidates, key=lambda item: item[0].name)
            data_origin[label] = {
                "manifest": path.name,
                "input_kind": manifest["input_kind"],
                "slurm_job_id": manifest["slurm_job_id"],
                "source_sha256": manifest["source_sha256"],
            }
    matched = {}
    _matching_stage(manifests, matched, checks, "03_fit_models", config_hash, [index_path])
    _matching_stage(manifests, matched, checks, "04_conformal", config_hash, [calibration_path])
    plan_path = outputs / "evaluation/plan.json"
    plan_binding = None
    if plan_path.exists():
        plan = read_json(plan_path)
        checks.check(plan["config_sha256"] == config_hash, "frozen plan config binding")
        checks.digest(index_path, plan["fit_index_sha256"], "frozen plan fit-index binding")
        checks.digest(
            calibration_path, plan["calibration_sha256"], "frozen plan calibration binding"
        )
        checks.check(
            plan["test_partition_sha256"] == partition_metadata["partitions_sha256"]["test"],
            "frozen plan test-partition metadata binding",
        )
        for relative, expected in plan["analysis_source_sha256"].items():
            checks.digest(root / relative, expected, f"frozen analysis source: {relative}")
        plan_binding = {
            "status": plan["status"],
            "slurm_job_id": plan.get("slurm_job_id"),
            "analysis_files": len(plan["analysis_source_sha256"]),
        }
        evaluation_artifacts = [plan_path] + [
            outputs / "metrics" / name
            for name in (
                "performance.json",
                "conformal_coverage.json",
                "rq2_feature_set_comparison.json",
            )
        ]
        _matching_stage(
            manifests, matched, checks, "05_evaluate", config_hash, evaluation_artifacts
        )
    elif require_exports:
        checks.check(False, "frozen evaluation plan exists")

    export_candidates = [
        (p, m)
        for p, m in manifests
        if m.get("stage") == "06_export_results"
        and m.get("config_sha256") == config_hash
        and m.get("status") == "complete"
    ]
    if export_candidates:
        matched["06_export_results"] = max(export_candidates, key=lambda item: item[0].name)
    elif require_exports:
        checks.check(False, "successful stage06 export manifest exists")

    source_digests = {
        record["source_sha256"] for record in (index["models"] | index["sensitivities"]).values()
    } | {manifest["source_sha256"] for _, manifest in matched.values()}
    source_digests |= {row["source_sha256"] for row in data_origin.values()}
    archive_members = {}
    config_relative = str(config_path.resolve().relative_to(root.resolve()))
    for digest in sorted(source_digests):
        path = outputs / "manifests" / f"source_{digest}.tar.gz"
        members = archive_hashes(path, digest, checks)
        archive_members[digest] = members
        checks.check(members.get(config_relative) == config_hash, f"archived config: {digest}")
        checks.check(members.get("uv.lock") == lock_hash, f"archived lock: {digest}")
    for key, record in (index["models"] | index["sensitivities"]).items():
        archived = archive_members.get(record["source_sha256"], {})
        for path, expected in index["fit_source_sha256"].items():
            relative = str(Path(path).relative_to(root))
            checks.check(
                archived.get(relative) == expected, f"fitted source snapshot: {key}:{relative}"
            )

    export_source_files = []
    if "06_export_results" in matched:
        export_manifest = matched["06_export_results"][1]
        archived = archive_members.get(export_manifest["source_sha256"], {})
        export_source_files = [
            *sorted((root / "src/brca/viz").glob("*.py")),
            root / "src/brca/reporting.py",
            root / "scripts/06_export_results.py",
        ]
        for path in export_source_files:
            relative = str(path.relative_to(root))
            checks.check(
                archived.get(relative) == sha256_file(path),
                f"current export source matches producing snapshot: {relative}",
            )

    provenance = {}
    for stage, (path, manifest) in matched.items():
        for field in (
            "git_sha",
            "git_dirty",
            "config_sha256",
            "seed",
            "package_versions",
            "slurm_job_id",
            "wall_clock_seconds",
            "source_sha256",
            "source_archive",
        ):
            checks.check(
                field in manifest and manifest[field] is not None,
                f"manifest required field: {stage}:{field}",
            )
        checks.check(manifest["seed"] == cfg.seed, f"manifest seed: {stage}")
        result = subprocess.run(
            ["git", "-C", str(root), "cat-file", "-e", manifest["git_sha"] + "^{commit}"],
            capture_output=True,
            check=False,
        )
        checks.check(result.returncode == 0, f"manifest Git commit exists: {stage}")
        mismatches = [
            f"{name}=={version}"
            for name, version in manifest["package_versions"].items()
            if name != "python" and version not in locked.get(normalize_package(name), set())
        ]
        checks.check(not mismatches, f"manifest packages match lock: {stage}:{mismatches}")
        for artifact, expected in manifest["artifacts_sha256"].items():
            checks.digest(artifact, expected, f"{stage} output digest: {Path(artifact).name}")
        for artifact, expected in manifest.get("input_artifacts_sha256", {}).items():
            checks.digest(artifact, expected, f"{stage} input digest: {Path(artifact).name}")
        provenance[stage] = {
            key: manifest[key]
            for key in (
                "git_sha",
                "git_dirty",
                "allow_dirty",
                "slurm_job_id",
                "source_sha256",
                "wall_clock_seconds",
                "status",
            )
        }
        provenance[stage].update(
            manifest=path.name,
            recorded_packages=len(manifest["package_versions"]) - 1,
            output_artifacts=len(manifest["artifacts_sha256"]),
            input_artifacts=len(manifest.get("input_artifacts_sha256", {})),
        )

    dirty = any(manifest["git_dirty"] for _, manifest in matched.values())
    return {
        "schema_version": 1,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "audit_slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "audit_source_sha256": sha256_file(Path(__file__)),
        "config_sha256": config_hash,
        "status": "failed" if checks.failures else "integrity_checks_passed",
        "checks": checks.count,
        "errors": checks.failures,
        "models": len(index["models"]),
        "sensitivities": len(index["sensitivities"]),
        "hashed_registered_artifacts": len(index["artifacts_sha256"]),
        "source_archives_checked": len(source_digests),
        "export_source_files_checked": len(export_source_files),
        "calibration_states": len(calibration["states"]),
        "calibration_methods": dict(
            Counter(r["registered_method"] for r in calibration["states"].values())
        ),
        "calibration_diagnostic_statuses": dict(
            Counter(
                r["calibration"]["diagnostics"].get("status", "not_applicable_primary_envelope")
                for r in calibration["states"].values()
            )
        ),
        "producing_manifests": provenance,
        "data_origin_metadata": data_origin,
        "frozen_evaluation_plan": plan_binding,
        "clean_commit_reproduction": False,
        "dirty_source_snapshots": dirty,
        "numerical_rerun_performed": False,
        "scope": (
            "Patient partitions, evaluation caches and posterior/model objects were not decoded. "
            "Registered artifact bytes were hashed; partition digests were compared in metadata. "
            "This checks integrity, not numerical reproducibility. Uncommitted snapshots are "
            "explicitly recorded and do not satisfy a clean-commit reproduction requirement."
        ),
        "portability_note": (
            "Manifests and fit indices contain absolute artifact paths. Restore the original "
            "directory layout or explicitly migrate/revalidate paths when moving the analysis."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-exports", action="store_true")
    args = parser.parse_args()
    require_compute_node("result provenance audit")
    if not os.environ.get("SLURM_JOB_ID"):
        raise LoginNodeError("The provenance audit must run inside a SLURM allocation.")
    config_path = args.config.resolve()
    cfg = load_config(config_path)
    output = args.output or PROJECT_ROOT / cfg.paths.outputs / "reports/provenance_audit.json"
    report = audit_results(PROJECT_ROOT, config_path, require_exports=args.require_exports)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    raise SystemExit(bool(report["errors"]))


if __name__ == "__main__":
    main()
