"""Aggregate immutable-source audit before defining the real analysis cohort."""

from __future__ import annotations

import argparse

from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT
from brca.data.load import load_raw_clinical
from brca.data.prepare import raw_audit
from brca.manifest import sha256_file
from brca.pipeline import StageRun, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    require_compute_node("data audit")
    run = StageRun("00b_data_audit", args.config, allow_dirty=args.allow_dirty)
    raw_path = PROJECT_ROOT / run.config.paths.raw / run.config.data_source.raw_subdir
    report = raw_audit(load_raw_clinical(raw_path), run.config)
    target = run.outputs / "metrics" / "data_audit.json"
    write_json(target, report)
    run.finish(
        [target],
        input_kind="real_metabric",
        raw_provenance_sha256=sha256_file(raw_path / "PROVENANCE.json"),
    )
    print(f"Aggregate data audit written: {target}")


if __name__ == "__main__":
    main()
