"""Stage 00 -- download the public METABRIC clinical tables from cBioPortal.

Network I/O only, no compute: safe on the login node. Writes an immutable raw
snapshot to ``data/raw/<provider>_<study>/`` and records its checksums in
``data/raw/CHECKSUMS.sha256``. See docs/metabric-data-dictionary.md.
"""

from __future__ import annotations

import argparse

from brca.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT, load_config
from brca.data.cbioportal import fetch_study, save_raw, to_wide, update_checksum_file


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--force", action="store_true", help="replace an existing raw snapshot")
    args = parser.parse_args()
    cfg = load_config(args.config)

    src = cfg.data_source
    raw_root = PROJECT_ROOT / cfg.paths.raw
    out_dir = raw_root / src.raw_subdir

    payloads = fetch_study(src.study_id)
    # Fail before writing anything if the response is not what the pipeline expects.
    patients = to_wide(payloads["patient"], "patientId")
    samples = to_wide(payloads["sample"], "sampleId")
    missing = sorted(set(src.required_patient_attributes) - set(patients.columns))
    missing += sorted(set(src.required_sample_attributes) - set(samples.columns))
    if missing:
        raise SystemExit(f"API response lacks required attributes: {missing}")

    checksums = save_raw(payloads, out_dir, force=args.force)
    update_checksum_file(raw_root / "CHECKSUMS.sha256", checksums, src.raw_subdir)
    rel = out_dir.relative_to(PROJECT_ROOT)
    print(f"saved {len(patients)} patients, {len(samples)} samples -> {rel}")


if __name__ == "__main__":
    main()
