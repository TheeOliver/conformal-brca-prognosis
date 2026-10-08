"""Explore the frozen training set and export aggregate figures."""

from __future__ import annotations

import argparse
import json

from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH
from brca.data.eda import training_eda
from brca.data.partitions import load_partition
from brca.manifest import sha256_file
from brca.pipeline import StageRun, write_json
from brca.viz.eda import export_training_eda


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    require_compute_node("training EDA")
    run = StageRun("02b_training_eda", args.config, allow_dirty=args.allow_dirty)
    cohort = load_partition(run.config, "train")
    metrics = run.outputs / "metrics" / "training_eda.json"
    write_json(metrics, training_eda(cohort, run.config))
    # Rendering reads only the saved aggregate; it never recomputes cohort metrics.
    figures = export_training_eda(
        json.loads(metrics.read_text()), run.outputs / "figures", run.config
    )
    run.finish(
        [metrics, *figures],
        input_kind="real_metabric",
        split_sha256=sha256_file(run.processed / "splits.json"),
        n_training=len(cohort),
    )
    print(f"Training EDA: {len(cohort)} patients, {len(figures)} figure files; {metrics}")


if __name__ == "__main__":
    main()
