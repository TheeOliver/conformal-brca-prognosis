"""Scaffold tests: the invariants the rules and hooks depend on.

These run in the Stop-hook gate, so they must stay fast and must not need real
data. Extend this suite as real analysis code lands -- especially the conformal
coverage assertions required by .claude/rules/conformal.md.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from brca.config import DEFAULT_CONFIG_PATH, load_config
from brca.data.load import load_cohort
from brca.data.schema import PAM50_LEVELS, SchemaError, validate

FIXTURE = Path(__file__).parent / "fixtures" / "synthetic_metabric.csv"


@pytest.fixture(scope="module")
def cfg():
    return load_config()


def test_config_file_exists():
    assert DEFAULT_CONFIG_PATH.exists(), "config/default.yaml is the single source of truth"


def test_split_fractions_sum_to_one(cfg):
    total = cfg.split.train + cfg.split.calibration + cfg.split.test
    assert total == pytest.approx(1.0), f"split fractions must sum to 1, got {total}"


def test_calibration_set_is_non_trivial(cfg):
    # Conformal quantiles need enough calibration points for the smallest alpha:
    # ceil((n+1)(1-alpha)) must not exceed n. See .claude/rules/conformal.md.
    assert cfg.split.calibration > 0.0


def test_alphas_are_miscoverage(cfg):
    alphas = cfg.conformal.alphas
    assert alphas, "at least one alpha required"
    assert all(0.0 < a < 0.5 for a in alphas), f"alpha is MISCOVERAGE, not coverage; got {alphas}"


def test_time_grid_percentiles_are_interior(cfg):
    pct = cfg.evaluation.time_grid_percentiles
    assert min(pct) > 0 and max(pct) < 100, (
        "time grid must be strictly interior or IPCW weights explode"
    )


def test_mcmc_gate_thresholds_present(cfg):
    d = cfg.mcmc.diagnostics
    assert d.max_rhat <= 1.01
    assert d.min_ess_bulk >= 400
    assert d.min_ess_tail >= 400
    assert d.max_divergences == 0


def test_synthetic_fixture_exists_and_validates():
    assert FIXTURE.exists(), "run: uv run python tests/fixtures/make_synthetic_metabric.py"
    df = load_cohort(FIXTURE)
    assert len(df) > 0


def test_fixture_retains_censored_rows():
    df = load_cohort(FIXTURE)
    assert not df["event"].all(), "censored rows must be retained, never dropped"
    assert df["event"].any(), "some events must be observed"


def test_fixture_covers_every_pam50_level():
    df = load_cohort(FIXTURE)
    assert set(df["pam50_subtype"]) == set(PAM50_LEVELS)


def test_fixture_is_deterministic():
    from tests.fixtures.make_synthetic_metabric import make_cohort

    assert make_cohort().equals(make_cohort())


def test_schema_rejects_integer_event_column():
    df = load_cohort(FIXTURE).assign(event=lambda d: d["event"].astype(int))
    with pytest.raises(SchemaError, match="must be bool"):
        validate(df)


def test_schema_rejects_non_positive_time():
    df = load_cohort(FIXTURE)
    df.loc[df.index[0], "survival_months"] = 0.0
    with pytest.raises(SchemaError, match="strictly positive"):
        validate(df)


def test_schema_rejects_duplicate_patients():
    df = load_cohort(FIXTURE)
    df.loc[df.index[1], "patient_id"] = df.loc[df.index[0], "patient_id"]
    with pytest.raises(SchemaError, match="duplicate"):
        validate(df)
