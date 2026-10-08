"""Independent aggregate/source audit; never deserialize patient-level artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tarfile
from collections import Counter
from pathlib import Path

import yaml


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", "--root", dest="root", type=Path, required=True)
    parser.add_argument("--job", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this independent audit through SLURM")
    root = args.root.resolve()
    outputs = root / "outputs"
    directory = outputs / "uncertainty_reliability"
    checks, failures = [], []

    def check(name, condition):
        checks.append(name)
        if not condition:
            failures.append(name)

    def hash_check(name, path, expected):
        check(name, Path(path).is_file() and digest(path) == expected)

    def near(a, b):
        if a is None or b is None:
            return a is b
        return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)

    manifests = [(p, read(p)) for p in (outputs / "manifests").glob("*.json")]
    runs = [
        (p, m)
        for p, m in manifests
        if m.get("stage") == "05_uncertainty_reliability" and str(m.get("slurm_job_id")) == args.job
    ]
    check("one completed experiment manifest", len(runs) == 1)
    if len(runs) != 1:
        raise RuntimeError(f"Expected one experiment manifest for job {args.job}")
    manifest_path, manifest = runs[0]
    check("experiment status", manifest["status"] == "exploratory_complete")
    check("experiment reportable", manifest.get("reportable") is True)
    check(
        "historical reconstruction ran",
        manifest.get("historical_point_statistics_verified") is True,
    )
    check("manifest interval cells", manifest["n_interval_cells"] == 96)
    check("manifest test count", manifest["n_test"] == 495)
    for path, sha in manifest["artifacts_sha256"].items():
        hash_check(f"output {Path(path).name}", path, sha)
    hash_check("source archive bytes", manifest["source_archive"], manifest["source_sha256"])

    plan = read(directory / "plan.json")
    prepared = read(directory / "preparation.json")
    binding = read(directory / "preparation_binding.json")
    coverage = read(directory / "coverage.json")
    method_pairs = read(directory / "paired_method_comparisons.json")
    feature_pairs = read(directory / "paired_feature_comparisons.json")
    performance = read(directory / "ridge_performance.json")
    index_path = outputs / "models/fit_index.json"
    index = read(index_path)
    config_path = Path(manifest["config_path"])
    config = yaml.safe_load(config_path.read_text())
    settings_path = root / "config/uncertainty_reliability.yaml"
    settings = yaml.safe_load(settings_path.read_text())
    processed = root / config["paths"]["processed"]

    check("post hoc status explicit", plan["status"] == "post_hoc_exploratory_test_reuse")
    check("settings contents frozen", plan["settings"] == settings)
    check("1000 planned ridge refits", settings["ridge_bootstrap_replicates"] == 1000)
    check("positive ridge penalty", settings["ridge_alpha"] > 0)
    check("config binding agrees", manifest["config_sha256"] == plan["config_sha256"])
    hash_check("configuration", config_path, plan["config_sha256"])
    hash_check("exploratory settings", settings_path, plan["settings_sha256"])
    hash_check("lock", root / "uv.lock", plan["lock_sha256"])
    hash_check("fit index", index_path, plan["fit_index_sha256"])
    hash_check("frozen split", processed / "splits.json", plan["split_sha256"])
    hash_check(
        "frozen test cache bytes only",
        outputs / "evaluation/test_cache.pkl",
        plan["test_cache_sha256"],
    )
    hash_check("original plan", outputs / "evaluation/plan.json", plan["original_plan_sha256"])
    for name, sha in plan["partition_sha256"].items():
        hash_check(f"partition bytes only {name}", processed / "partitions" / f"{name}.csv", sha)
    for path, sha in plan["reference_artifacts_sha256"].items():
        hash_check(f"reference {path}", path, sha)
    for path, sha in plan["analysis_source_sha256"].items():
        hash_check(f"current analysis source {path}", root / path, sha)
    with tarfile.open(manifest["source_archive"], "r:gz") as archive:
        for path, sha in plan["analysis_source_sha256"].items():
            entry = archive.extractfile(path)
            check(
                f"archived analysis source {path}",
                entry is not None and hashlib.sha256(entry.read()).hexdigest() == sha,
            )
        archived_config = archive.extractfile(str(config_path.relative_to(root)))
        check(
            "archived default config",
            hashlib.sha256(archived_config.read()).hexdigest() == plan["config_sha256"],
        )
        archived_settings = archive.extractfile("config/uncertainty_reliability.yaml")
        check(
            "archived exploratory settings",
            hashlib.sha256(archived_settings.read()).hexdigest() == plan["settings_sha256"],
        )
    hash_check("preparation plan binding", directory / "plan.json", binding["plan_sha256"])
    hash_check(
        "preparation state binding", directory / "preparation.json", binding["preparation_sha256"]
    )
    for path, sha in prepared["artifacts_sha256"].items():
        hash_check(f"prepared artifact {path}", path, sha)
    check("training size", prepared["n_train"] == 989)
    check("calibration size", prepared["n_cal"] == 495)

    original_paths = [
        outputs / "evaluation/plan.json",
        outputs / "metrics/conformal_coverage.json",
        outputs / "qin_exploratory/plan.json",
        outputs / "qin_exploratory/coverage.json",
    ]
    baseline_evidence = {}
    historical_sources_checked = set()
    for path in original_paths:
        earlier = [
            (p, m, m["artifacts_sha256"][str(path)])
            for p, m in manifests
            if m.get("stage") in {"05_evaluate", "05_qin_exploratory"}
            and str(path) in m.get("artifacts_sha256", {})
            and m.get("status") in {"complete", "post_hoc_exploratory_complete"}
        ]
        check(f"original manifest exists {path}", bool(earlier))
        matching = [(p, m, sha) for p, m, sha in earlier if digest(path) == sha]
        check(f"original manifest hash preserved {path}", bool(matching))
        if matching:
            original_manifest_path, original_manifest, sha = matching[0]
            baseline_evidence[str(path)] = {
                "sha256": sha,
                "original_manifest": str(original_manifest_path),
                "original_manifest_sha256": digest(original_manifest_path),
                "original_slurm_job_id": original_manifest["slurm_job_id"],
                "original_source_archive": original_manifest["source_archive"],
                "original_source_archive_sha256": original_manifest["source_sha256"],
            }
            if original_manifest_path not in historical_sources_checked:
                historical_sources_checked.add(original_manifest_path)
                hash_check(
                    f"historical source archive {original_manifest_path.name}",
                    original_manifest["source_archive"],
                    original_manifest["source_sha256"],
                )
                old_plan_path = outputs / (
                    "evaluation/plan.json"
                    if original_manifest["stage"] == "05_evaluate"
                    else "qin_exploratory/plan.json"
                )
                old_plan = read(old_plan_path)
                with tarfile.open(original_manifest["source_archive"], "r:gz") as archive:
                    for source_path, source_sha in old_plan["analysis_source_sha256"].items():
                        entry = archive.extractfile(source_path)
                        check(
                            f"historical source {original_manifest_path.name}/{source_path}",
                            entry is not None
                            and hashlib.sha256(entry.read()).hexdigest() == source_sha,
                        )

    def interval_key(model, alpha, horizon, method):
        return f"{model}__alpha{alpha:g}__horizon{horizon:g}__{method}"

    lower_method, ridge_method = "conservative_lower_bound", "qin_bootstrap_ridge_cox"
    expected_lower = {
        interval_key(key, alpha, horizon, lower_method)
        for key in index["models"]
        for alpha in plan["alphas"]
        for horizon in plan["horizons_months"]
    }
    cox_keys = {key for key, value in index["models"].items() if value["model"] == "cox"}
    expected_ridge = {
        interval_key(key, alpha, horizon, ridge_method)
        for key in cox_keys
        for alpha in plan["alphas"]
        for horizon in plan["horizons_months"]
    }
    check("72 planned lower cells", len(expected_lower) == 72)
    check("24 planned ridge cells", len(expected_ridge) == 24)
    check("all lower states present", set(prepared["lower_states"]) == expected_lower)
    check("all coverage cells present", set(coverage) == expected_lower | expected_ridge)
    check("all paired method cells present", set(method_pairs) == set(coverage))
    check("four ridge model arms", set(prepared["ridge"]) == cox_keys and len(cox_keys) == 4)
    ridge_counts = {}
    for key, entry in prepared["ridge"].items():
        diagnostics = entry["diagnostics"]
        ridge_counts[key] = diagnostics["bootstrap_successful"]
        for name in ("bootstrap_requested", "bootstrap_attempted", "bootstrap_successful"):
            check(f"{key} {name}", diagnostics[name] == 1000)
        check(f"{key} zero failed refits", diagnostics["bootstrap_failed"] == 0)
        check(f"{key} every pivot present", len(entry["pivots"]) == 1000)
        check(f"{key} valid pivots", all(math.isfinite(v) and 0 <= v <= 1 for v in entry["pivots"]))
        check(
            f"{key} matched original/refit parameters",
            diagnostics["original_bootstrap_parameter_match"] is True,
        )
        check(
            f"{key} ridge penalty frozen",
            diagnostics["working_model_parameters"]["alpha"] == settings["ridge_alpha"],
        )

    original_coverage = read(outputs / "metrics/conformal_coverage.json")
    original_qin = read(outputs / "qin_exploratory/coverage.json")
    minimum_group = config["eda"]["minimum_group_size"]
    suppressed_count = 0
    summary_rows = []
    comparable = [
        "observable_coverage_lower",
        "observable_coverage_upper",
        "width_median_months",
        "width_q25_months",
        "width_q75_months",
        "full_support_fraction",
    ]
    for key, entry in coverage.items():
        marginal = entry["marginal"]
        is_lower = entry["method"] == lower_method
        check(f"{key} n_test", marginal["n_test"] == 495)
        check(f"{key} n_cal", marginal["n_cal"] == (495 if is_lower else 1000))
        check(
            f"{key} coverage ordered",
            0
            <= marginal["observable_coverage_lower"]
            <= marginal["observable_coverage_upper"]
            <= 1,
        )
        check(
            f"{key} width ordered",
            0
            <= marginal["width_q25_months"]
            <= marginal["width_median_months"]
            <= marginal["width_q75_months"]
            <= marginal["horizon_months"],
        )
        for column, groups in entry["subgroups"].items():
            check(
                f"{key} subgroup counts {column}",
                sum(value["n_test"] for value in groups.values()) == 495,
            )
            for level, value in groups.items():
                if value["n_test"] < minimum_group:
                    suppressed_count += 1
                    check(
                        f"{key} suppression {column}/{level}",
                        value == {"n_test": value["n_test"], "status": "suppressed_small_group"},
                    )
                else:
                    check(
                        f"{key} unsuppressed {column}/{level}",
                        value.get("status") != "suppressed_small_group"
                        and "width_median_months" in value,
                    )
        if is_lower:
            old_key = key.removesuffix(lower_method) + "conservative"
            reference = original_coverage[old_key]["marginal"]
        else:
            old_key = key.removesuffix(ridge_method) + "qin_bootstrap_cox"
            reference = original_qin[old_key]["marginal"]
        for metric in comparable:
            check(
                f"{key} paired point {metric}",
                near(
                    method_pairs[key]["differences"][metric]["estimate"],
                    marginal[metric] - reference[metric],
                ),
            )
        if marginal["alpha"] == settings["report_alpha"]:
            summary_rows.append(
                {
                    "key": key,
                    **{
                        name: marginal[name]
                        for name in [
                            "alpha",
                            "horizon_months",
                            "observable_coverage_lower",
                            "observable_coverage_upper",
                            "ipcw_coverage",
                            "width_median_months",
                        ]
                    },
                    "width_difference": method_pairs[key]["differences"]["width_median_months"],
                }
            )
    clinical_keys = {key for key in coverage if "__clinical__" in key}
    check(
        "48 paired feature cells", set(feature_pairs) == clinical_keys and len(feature_pairs) == 48
    )
    for key, entry in feature_pairs.items():
        other = key.replace("__clinical__", "__clinical_molecular__")
        for metric in comparable:
            check(
                f"{key} paired feature point {metric}",
                near(
                    entry["differences"][metric]["estimate"],
                    coverage[other]["marginal"][metric] - coverage[key]["marginal"][metric],
                ),
            )
    check("both endpoints evaluated", set(performance) == {"disease_specific", "overall"})
    for endpoint, entry in performance.items():
        check(f"{endpoint} performance reportable", entry["reportable"] is True)
        check(
            f"{endpoint} both ridge feature arms",
            set(entry["models"]) == {"cox_ridge__clinical", "cox_ridge__clinical_molecular"},
        )

    def no_gate_failure(value):
        if isinstance(value, dict):
            return value.get("status") != "bootstrap_gate_failed" and all(
                no_gate_failure(item) for item in value.values()
            )
        if isinstance(value, list):
            return all(no_gate_failure(item) for item in value)
        return True

    check(
        "no hidden bootstrap gate failures",
        no_gate_failure([coverage, method_pairs, feature_pairs, performance]),
    )
    result = {
        "status": "passed" if not failures else "failed",
        "audit_slurm_job_id": os.environ["SLURM_JOB_ID"],
        "analysis_slurm_job_id": args.job,
        "analysis_manifest": str(manifest_path),
        "analysis_manifest_sha256": digest(manifest_path),
        "audit_source": str(Path(__file__).resolve()),
        "audit_source_sha256": digest(__file__),
        "baseline_evidence": baseline_evidence,
        "checks_passed": len(checks) - len(failures),
        "checks_total": len(checks),
        "failures": failures,
        "coverage_cells": dict(Counter(row["method"] for row in coverage.values())),
        "ridge_successful_refits": ridge_counts,
        "suppressed_subgroup_cells": suppressed_count,
        "ninety_percent_summary": summary_rows,
        "scope": "source_artifact_and_aggregate_consistency_not_numerical_reproduction",
        "patient_artifacts": "hashed_only_never_deserialized",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "checks_passed",
                    "checks_total",
                    "failures",
                    "coverage_cells",
                    "ridge_successful_refits",
                )
            },
            indent=2,
        )
    )
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
