"""The login-node guard is the one compute safeguard every runner is subject to."""

from __future__ import annotations

import pytest

from brca.compute_guard import LoginNodeError, on_login_node, require_compute_node

LOGIN = "openlab-slurm.hpc.local"
COMPUTE = "gpgpu01"


def test_login_node_detected_outside_allocation():
    assert on_login_node(LOGIN, env={})


def test_compute_node_is_not_login_node():
    assert not on_login_node(COMPUTE, env={})


def test_inside_slurm_allocation_is_allowed_even_on_login_host():
    assert not on_login_node(LOGIN, env={"SLURM_JOB_ID": "12345"})


def test_require_raises_on_login_node_and_names_sbatch():
    with pytest.raises(LoginNodeError, match="sbatch"):
        require_compute_node("stage 03", hostname=LOGIN, env={})


def test_require_passes_inside_allocation():
    require_compute_node("stage 03", hostname=LOGIN, env={"SLURM_JOB_ID": "1"})


def test_explicit_override_is_honoured():
    require_compute_node("stage 03", hostname=LOGIN, env={"BRCA_ALLOW_LOGIN_NODE": "1"})
