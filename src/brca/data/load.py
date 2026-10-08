"""The single entry point for reading cohort data.

Never bypass this with a bare ``pd.read_csv`` -- validation is the point.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd

from brca.data.cbioportal import RAW_FILES, to_wide
from brca.data.schema import validate
from brca.manifest import sha256_file


def load_cohort(path: Path | str) -> pd.DataFrame:
    """Read a cohort table and validate it against the schema contract."""
    # IDs are identifiers, including when they look numeric (e.g. leading zeros).
    return validate(pd.read_csv(path, dtype={"patient_id": str}))


def load_raw_clinical(path: Path | str) -> pd.DataFrame:
    """Validate the immutable source checksums and join one sample per patient.

    This returns source column names for stage 01 and aggregate source auditing.
    Multiple samples per patient require an explicit cohort rule, never a silent pick.
    """
    root = Path(path)
    provenance = json.loads((root / "PROVENANCE.json").read_text())
    if not {RAW_FILES["patient"], RAW_FILES["sample"]} <= set(provenance["files_sha256"]):
        raise ValueError("Raw snapshot provenance lacks required clinical checksums.")
    for filename, expected in provenance["files_sha256"].items():
        if Path(filename).name != filename or sha256_file(root / filename) != expected:
            raise ValueError("Raw snapshot checksum mismatch; restore the original snapshot.")
    payloads = {}
    for key in ("patient", "sample"):
        with gzip.open(root / RAW_FILES[key], "rt") as handle:
            payloads[key] = json.load(handle)
    patients = to_wide(payloads["patient"], "patientId")
    samples = pd.DataFrame(payloads["sample"])[["patientId", "sampleId"]].drop_duplicates()
    if samples["patientId"].duplicated().any() or samples["sampleId"].duplicated().any():
        raise ValueError("Raw clinical data does not have exactly one sample per patient.")
    sample_values = to_wide(payloads["sample"], "patientId")
    overlap = set(patients.columns) & set(sample_values.columns)
    if overlap:
        raise ValueError("Patient and sample tables have ambiguous overlapping attributes.")
    result = patients.join(sample_values, how="left", validate="one_to_one")
    result.index.name = "patient_id"
    return result
