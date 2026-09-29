"""Stage 01 -- load the raw cohort, validate the schema, write the processed table.

Stage stub. Implement against the rules that scope this path -- start with
.claude/rules/splits-and-leakage.md and .claude/rules/reproducibility.md.
"""

from __future__ import annotations

import argparse

from brca.config import DEFAULT_CONFIG_PATH, load_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    args = parser.parse_args()
    load_config(args.config)
    raise NotImplementedError(
        "scripts/01_prepare_data.py is a stub; see docs/experimental-design.md"
    )


if __name__ == "__main__":
    main()
