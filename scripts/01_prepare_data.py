"""Stage 01 -- load the raw cohort, validate the schema, write the processed table.

Maps the checksummed public clinical snapshot into the registered DSS cohort,
retaining OS labels for sensitivity analysis. No imputation is fitted here.
"""

from __future__ import annotations

import argparse

from brca.compute_guard import require_compute_node
from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT
from brca.data.load import load_raw_clinical
from brca.data.partitions import archive_synthetic_split
from brca.data.prepare import prepare_cohort
from brca.manifest import sha256_file
from brca.pipeline import StageRun, write_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument(
        "--archive-synthetic",
        action="store_true",
        help="preserve an existing synthetic split before real-data preparation",
    )
    args = parser.parse_args()
    require_compute_node("stage 01")
    run = StageRun("01_prepare_data", args.config, allow_dirty=args.allow_dirty)
    if args.archive_synthetic:
        archive_synthetic_split(run.processed, run.outputs)
    raw_path = PROJECT_ROOT / run.config.paths.raw / run.config.data_source.raw_subdir
    cohort, report = prepare_cohort(load_raw_clinical(raw_path), run.config)
    run.processed.mkdir(parents=True, exist_ok=True)
    target = run.processed / "cohort.csv"
    content = cohort.to_csv(index=False, lineterminator="\n").encode()
    if target.exists() and target.read_bytes() != content:
        raise ValueError(
            "Existing prepared cohort differs; archive the experiment before replacing it."
        )
    if not target.exists():
        target.write_bytes(content)
    report_path = run.outputs / "metrics" / "cohort_preparation.json"
    write_json(report_path, report)
    run.finish(
        [target, report_path],
        input_kind="real_metabric",
        raw_provenance_sha256=sha256_file(raw_path / "PROVENANCE.json"),
    )
    print(f"Prepared cohort: {len(cohort)} patients; exclusion counts recorded in {report_path}")


if __name__ == "__main__":
    main()
