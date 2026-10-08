"""Stage 04: calibrate frozen training models using calibration patients only."""

import argparse

from brca.calibration import run_calibration
from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH
from brca.pipeline import StageRun


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    require_compute_node("stage 04 (conformal)")
    run = StageRun("04_conformal", args.config, allow_dirty=args.allow_dirty)
    print(f"Calibration: {run_calibration(run)}")


if __name__ == "__main__":
    main()
