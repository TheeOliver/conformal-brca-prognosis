"""Session-wide test setup shared by every agent and human running pytest."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from brca.compute_guard import OVERRIDE_ENV, on_login_node

FIXTURE = Path(__file__).parent / "fixtures" / "synthetic_metabric.csv"


def pytest_sessionstart(session: pytest.Session) -> None:
    # The fixture is generated, not committed: no .csv ever enters git, which
    # keeps the data rule absolute and simple to enforce.
    if not FIXTURE.exists():
        from tests.fixtures.make_synthetic_metabric import main

        main()


@pytest.hookimpl(trylast=True)  # after `-m` deselection has happened
def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    slow = [item for item in items if item.get_closest_marker("slow")]
    if slow and on_login_node() and os.environ.get(OVERRIDE_ENV) != "1":
        pytest.exit(
            f"{len(slow)} slow test(s) selected on the shared login node. They sample real "
            "chains -- run them under sbatch (docs/slurm-guide.md), or use `make test`.",
            returncode=3,
        )
