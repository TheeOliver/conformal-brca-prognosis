"""The single entry point for reading cohort data.

Never bypass this with a bare ``pd.read_csv`` -- validation is the point.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from brca.data.schema import validate


def load_cohort(path: Path | str) -> pd.DataFrame:
    """Read a cohort table and validate it against the schema contract."""
    # IDs are identifiers, including when they look numeric (e.g. leading zeros).
    return validate(pd.read_csv(path, dtype={"patient_id": str}))
