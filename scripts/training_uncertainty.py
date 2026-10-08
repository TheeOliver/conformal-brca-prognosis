"""Illustrate Bayesian parameter and individual-outcome uncertainty using training aggregates."""

from __future__ import annotations

import argparse
from pathlib import Path

from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH
from brca.data.partitions import load_partition
from brca.data.preprocessing import TrainingPreprocessor
from brca.fitting import load_fit, load_fit_index
from brca.manifest import sha256_file
from brca.pipeline import StageRun, write_json
from brca.uncertainty import hypothetical_profiles, hypothetical_uncertainty_report
from brca.viz.uncertainty import export_hypothetical_uncertainty


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    require_compute_node("training hypothetical-profile uncertainty")
    run = StageRun("training_uncertainty", args.config, allow_dirty=args.allow_dirty)
    index = load_fit_index(run)
    train = load_partition(run.config, "train")
    profiles, definitions = hypothetical_profiles(train, run.config)
    models, preprocessors = {}, {}
    for feature in ("clinical", "clinical_molecular"):
        key = f"{run.config.outcome.primary_endpoint}__weibull__{feature}"
        models[feature] = load_fit(index["models"][key])
        preprocessors[feature] = TrainingPreprocessor.load(Path(index["preprocessors"][feature]))
    report = hypothetical_uncertainty_report(
        profiles, definitions, models, preprocessors, run.config
    )
    report["fit_index_sha256"] = sha256_file(run.outputs / "models/fit_index.json")
    report["training_partition_sha256"] = index["training_sha256"]
    report["config_sha256"] = run.metadata["config_sha256"]
    report["n_training_for_aggregates"] = len(train)
    metric_path = run.outputs / "metrics/rq4_hypothetical_uncertainty.json"
    write_json(metric_path, report)
    figures = export_hypothetical_uncertainty(metric_path, run.outputs / "figures", run.config)
    run.finish(
        [metric_path, *figures],
        status="complete",
        n_hypothetical_profiles=len(profiles),
        training_only=True,
        calibration_accessed=False,
        test_accessed=False,
    )
    print(f"Hypothetical-profile uncertainty: {len(profiles)} profiles; {len(figures)} figures")


if __name__ == "__main__":
    main()
