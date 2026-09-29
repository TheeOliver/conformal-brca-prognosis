"""Stage 05 -- score on test ONCE and write outputs/metrics/. COMPUTE.

Stage stub. Implement against the rules that scope this path -- start with
.claude/rules/splits-and-leakage.md and .claude/rules/reproducibility.md.
"""

from __future__ import annotations

import argparse

from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH, load_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    args = parser.parse_args()
    require_compute_node("stage 05 (evaluate)")
    load_config(args.config)
    raise NotImplementedError("scripts/05_evaluate.py is a stub; see docs/experimental-design.md")


if __name__ == "__main__":
    main()
