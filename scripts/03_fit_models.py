"""Stage 03: fit the registered model matrix and sensitivities on training patients."""

import argparse

from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH
from brca.fitting import run_fitting
from brca.pipeline import StageRun


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    require_compute_node("stage 03 (fit)")
    run = StageRun("03_fit_models", args.config, allow_dirty=args.allow_dirty)
    print(f"Model index: {run_fitting(run)}")


if __name__ == "__main__":
    main()
