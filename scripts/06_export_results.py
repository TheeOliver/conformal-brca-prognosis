"""Stage 06: export thesis figures/tables from saved aggregate metrics only."""

import argparse
import json
from pathlib import Path

from brca.calibration import load_calibration
from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH
from brca.fitting import load_fit_index
from brca.manifest import sha256_file
from brca.pipeline import StageRun
from brca.reporting import write_research_report
from brca.viz.coverage import render_coverage
from brca.viz.eda import export_training_eda
from brca.viz.results import render_bayesian_workflow, render_performance


def validated_export_inputs(run):
    index = load_fit_index(run)
    calibration = load_calibration(run, index)
    calibration_path = run.outputs / "models/conformal_calibration.json"
    calibration_hash = sha256_file(calibration_path)
    calibration_verified = False
    for path in (run.outputs / "manifests").glob("04_conformal_*.json"):
        manifest = json.loads(path.read_text())
        if (
            manifest.get("status") == "complete"
            and manifest["config_sha256"] == run.metadata["config_sha256"]
            and manifest["artifacts_sha256"].get(str(calibration_path)) == calibration_hash
            and sha256_file(manifest["source_archive"]) == manifest["source_sha256"]
        ):
            calibration_verified = True
            break
    if not calibration_verified:
        raise ValueError("Calibration lacks a matching successful stage04 manifest.")
    plan = json.loads((run.outputs / "evaluation/plan.json").read_text())
    if (
        plan["config_sha256"] != run.metadata["config_sha256"]
        or plan["fit_index_sha256"] != calibration["fit_index_sha256"]
        or plan["calibration_sha256"]
        != sha256_file(run.outputs / "models/conformal_calibration.json")
    ):
        raise ValueError("Export inputs do not match the frozen evaluation plan.")
    paths = [
        run.outputs / "metrics" / name
        for name in (
            "performance.json",
            "conformal_coverage.json",
            "rq2_feature_set_comparison.json",
        )
    ]
    paths.append(run.outputs / "evaluation/plan.json")
    current_hashes = {str(path): sha256_file(path) for path in paths}
    for manifest_path in sorted(
        (run.outputs / "manifests").glob("05_evaluate_*.json"), reverse=True
    ):
        manifest = json.loads(manifest_path.read_text())
        if (
            manifest["config_sha256"] == run.metadata["config_sha256"]
            and all(manifest["artifacts_sha256"].get(p) == h for p, h in current_hashes.items())
            and sha256_file(manifest["source_archive"]) == manifest["source_sha256"]
        ):
            return index
    raise ValueError("Final metrics lack a matching source/config/artifact manifest.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    require_compute_node("stage 06 (export)")
    run = StageRun("06_export_results", args.config, allow_dirty=args.allow_dirty)
    metrics = run.outputs / "metrics"
    index = validated_export_inputs(run)
    artifacts = []
    artifacts += export_training_eda(
        json.loads((metrics / "training_eda.json").read_text()),
        run.outputs / "figures",
        run.config,
    )
    performance = json.loads((metrics / "performance.json").read_text())
    for endpoint, report in performance.items():
        artifacts += render_performance(report, run.outputs, f"rq1_rq2_rq3_{endpoint}", run.config)
    coverage = json.loads((metrics / "conformal_coverage.json").read_text())
    artifacts += render_coverage(coverage, run.outputs, run.config)
    for key, record in (index["models"] | index["sensitivities"]).items():
        if record["model"] == "weibull":
            report = json.loads((metrics / f"bayes_diagnostics_{key}.json").read_text())
            artifacts += render_bayesian_workflow(
                report["workflow"],
                report["diagnostics"],
                run.outputs,
                f"rq4_{key}",
                run.config,
            )
    artifacts += write_research_report(run.outputs, run.config)
    inputs = [
        metrics / name
        for name in (
            "performance.json",
            "conformal_coverage.json",
            "rq2_feature_set_comparison.json",
            "bayesian_sensitivity.json",
            "training_eda.json",
        )
    ] + [run.outputs / "models/fit_index.json"]
    inputs += [
        Path(record["metrics_path"])
        for record in (index["models"] | index["sensitivities"]).values()
    ]
    run.finish(
        artifacts,
        status="complete",
        input_artifacts_sha256={str(path): sha256_file(path) for path in inputs},
    )
    print(f"Exported {len(artifacts)} thesis artifacts.")


if __name__ == "__main__":
    main()
