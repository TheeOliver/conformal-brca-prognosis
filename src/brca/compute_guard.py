"""Refuse to run compute on the shared SLURM login node.

This guard lives in the code rather than in any agent's hook configuration, so
it holds for Claude, Codex and a human at a terminal alike. Agent hooks refuse
the same commands earlier, but only this layer is guaranteed to run.
See CLAUDE.md / AGENTS.md convention 1 and docs/slurm-guide.md.
"""

from __future__ import annotations

import os
import socket
from collections.abc import Mapping

LOGIN_HOST_MARKERS = ("slurm", "login", "head")
OVERRIDE_ENV = "BRCA_ALLOW_LOGIN_NODE"


class LoginNodeError(RuntimeError):
    """Raised when compute is attempted on the login node."""


def on_login_node(hostname: str | None = None, env: Mapping[str, str] | None = None) -> bool:
    """True when running on the login node outside any SLURM allocation."""
    env = os.environ if env is None else env
    if env.get("SLURM_JOB_ID"):
        return False
    host = (hostname or socket.gethostname()).lower()
    return any(marker in host for marker in LOGIN_HOST_MARKERS)


def require_compute_node(
    stage: str, hostname: str | None = None, env: Mapping[str, str] | None = None
) -> None:
    """Raise :class:`LoginNodeError` unless running inside a SLURM allocation.

    ``BRCA_ALLOW_LOGIN_NODE=1`` overrides, for a deliberate tiny smoke run only.
    """
    env = os.environ if env is None else env
    if env.get(OVERRIDE_ENV) == "1":
        return
    if on_login_node(hostname, env):
        raise LoginNodeError(
            f"{stage} is compute and must not run on the shared login node. "
            "Submit it instead: sbatch scripts/slurm/run_stage.sh <NN> "
            "(see docs/slurm-guide.md)."
        )
