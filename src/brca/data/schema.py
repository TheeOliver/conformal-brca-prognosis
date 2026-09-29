"""The cohort data contract.

Validation fails loudly. A dataset that does not meet the contract is a problem to
fix at the source, never something to coerce silently.
"""

from __future__ import annotations

import pandas as pd

REQUIRED_COLUMNS: tuple[str, ...] = (
    "patient_id",
    "survival_months",
    "event",
    "age_at_diagnosis",
    "tumor_size",
    "tumor_grade",
    "lymph_nodes_positive",
    "er_status",
    "pr_status",
    "her2_status",
    "pam50_subtype",
)

PAM50_LEVELS: tuple[str, ...] = ("LumA", "LumB", "Her2", "Basal", "Normal")


class SchemaError(ValueError):
    """Raised when a cohort does not meet the contract."""


def validate(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` unchanged, or raise :class:`SchemaError` explaining why not."""
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise SchemaError(f"missing required columns: {missing}")

    if df["event"].dtype != bool:
        raise SchemaError(
            "`event` must be bool, not 0/1 int -- see .claude/rules/survival-labels.md"
        )
    if not df["survival_months"].gt(0).all():
        raise SchemaError("`survival_months` must be strictly positive; do not clip")
    if df["patient_id"].duplicated().any():
        raise SchemaError("duplicate patient_id -- a leakage hazard across splits")
    if not df["event"].any():
        raise SchemaError("no observed events; every row is censored")
    if df["event"].all():
        raise SchemaError("no censored rows; censoring information appears to be lost")
    return df
