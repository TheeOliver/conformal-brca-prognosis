"""Stage 02 -- draw the ONE train/calibration/test split and freeze it.

Lightweight and safe on the login node. Use --synthetic for a development run
until stage 01 supplies data/processed/cohort.csv. Compatible reruns reuse the
frozen split; changed inputs never overwrite it.
"""

from __future__ import annotations

import argparse

from brca.config import DEFAULT_CONFIG_PATH
from brca.data.splits import run_split_stage


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--input", help="validated cohort CSV (default: processed/cohort.csv)")
    source.add_argument("--synthetic", action="store_true", help="use the synthetic test fixture")
    parser.add_argument(
        "--allow-dirty", action="store_true", help="permit an uncommitted development run"
    )
    args = parser.parse_args()
    try:
        split_path, manifest_path = run_split_stage(
            args.config,
            input_path=args.input,
            synthetic=args.synthetic,
            allow_dirty=args.allow_dirty,
        )
    except (ValueError, RuntimeError) as exc:
        parser.exit(1, f"Stage 02: {exc}\n")
    print(f"Frozen split: {split_path}\nRun manifest: {manifest_path}")


if __name__ == "__main__":
    main()
