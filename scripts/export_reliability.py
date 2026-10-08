"""Export audited, saved uncertainty-reliability aggregates on a SLURM worker."""

import argparse
import json

from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT
from brca.manifest import sha256_file
from brca.pipeline import StageRun, write_json
from brca.viz.reliability import render_reliability


def validated_reliability_inputs(run):
    """Read aggregate JSON only, bound to a successful frozen evaluation manifest."""
    directory = run.outputs / "uncertainty_reliability"
    paths = {
        name: directory / f"{name}.json"
        for name in (
            "coverage",
            "paired_method_comparisons",
            "paired_feature_comparisons",
            "ridge_performance",
            "preparation",
            "plan",
        )
    }
    inputs = {name: json.loads(path.read_text()) for name, path in paths.items()}
    plan = inputs["plan"]
    if plan["config_sha256"] != run.metadata["config_sha256"]:
        raise ValueError("Export config differs from the frozen reliability protocol")
    for relative, digest in plan["analysis_source_sha256"].items():
        if sha256_file(PROJECT_ROOT / relative) != digest:
            raise ValueError(f"Statistical source changed after the reliability run: {relative}")
    hashes = {str(path): sha256_file(path) for path in paths.values()}
    matched = None
    for candidate in sorted(
        (run.outputs / "manifests").glob("05_uncertainty_reliability_*.json"), reverse=True
    ):
        manifest = json.loads(candidate.read_text())
        if (
            manifest.get("status") == "exploratory_complete"
            and manifest.get("reportable") is True
            and manifest["config_sha256"] == plan["config_sha256"]
            and all(
                manifest["artifacts_sha256"].get(path) == digest for path, digest in hashes.items()
            )
            and sha256_file(manifest["source_archive"]) == manifest["source_sha256"]
        ):
            matched = candidate
            break
    if matched is None:
        raise ValueError("Reliability inputs lack a successful matching source/artifact manifest")
    references = {
        "primary_coverage": run.outputs / "metrics/conformal_coverage.json",
        "original_qin_coverage": run.outputs / "qin_exploratory/coverage.json",
        "primary_performance": run.outputs / "metrics/performance.json",
    }
    for name, path in references.items():
        if plan["reference_artifacts_sha256"].get(str(path)) != sha256_file(path):
            raise ValueError(f"Frozen reference changed: {name}")
        paths[name] = path
        inputs[name] = json.loads(path.read_text())
    return inputs, paths, matched


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    require_compute_node("reliability export")
    run = StageRun("06_export_reliability", args.config, allow_dirty=args.allow_dirty)
    inputs, paths, evaluation_manifest = validated_reliability_inputs(run)
    hashes = {str(path): sha256_file(path) for path in paths.values()}
    artifacts = render_reliability(inputs, run.outputs, run.config)
    if any(sha256_file(path) != digest for path, digest in hashes.items()):
        raise ValueError("An aggregate input changed during export")
    metadata = run.outputs / "uncertainty_reliability/export_metadata.json"
    write_json(
        metadata,
        {
            "status": "post_hoc_exploratory_test_reuse",
            "statistics_recomputed": False,
            "input_artifacts_sha256": hashes,
            "evaluation_manifest": str(evaluation_manifest),
            "evaluation_manifest_sha256": sha256_file(evaluation_manifest),
            "export_source_sha256": run.metadata["source_sha256"],
            "export_source_archive": run.metadata["source_archive"],
            "artifacts_sha256": {str(path): sha256_file(path) for path in artifacts},
            "figure_report_alpha": inputs["plan"]["settings"]["report_alpha"],
            "coverage_interval_meaning": "observable_bounds_are_not_confidence_intervals",
            "one_sided_target": "[L,tau] for min(T,tau), equivalently T>=L",
        },
    )
    run.finish(
        [*artifacts, metadata],
        status="complete",
        input_artifacts_sha256=hashes,
        evaluation_manifest_sha256=sha256_file(evaluation_manifest),
    )
    print(f"Exported {len(artifacts)} reliability artifacts plus provenance metadata.")


if __name__ == "__main__":
    main()
