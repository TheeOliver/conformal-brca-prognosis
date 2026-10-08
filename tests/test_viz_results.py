"""Publication exports preserve aggregate estimates and disclose unavailable CIs."""

from __future__ import annotations

from copy import deepcopy
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

from brca.viz.results import render_bayesian_workflow, render_performance


def figure_config():
    return SimpleNamespace(
        figures=SimpleNamespace(
            width_inches=6.3, height_inches=4, dpi=120, font_size=9, line_width=2
        ),
        bayesian=SimpleNamespace(predictive_interval_probability=0.90),
    )


def interval(estimate, available=True):
    return {
        "estimate": estimate,
        "ci_lower": estimate - 0.02 if available else None,
        "ci_upper": estimate + 0.02 if available else None,
        "status": "ok" if available else "bootstrap_gate_failed",
    }


def performance_fixture():
    def model(offset, available=True):
        return {
            "uno_c": interval(0.71 + offset, available),
            "ibs": interval(0.15 - offset),
            "mean_auc": interval(0.74 + offset),
            "brier": [
                {"time_months": t, **interval(v - offset)}
                for t, v in [(12, 0.08), (60, 0.15), (120, 0.22)]
            ],
            "auc": [
                {"time_months": t, **interval(v + offset, available)}
                for t, v in [(12, 0.76), (60, 0.74), (120, 0.72)]
            ],
        }

    return {
        "bootstrap": {"confidence_level": 0.95},
        "models": {"cox__clinical": model(0), "cox__clinical_molecular": model(0.01, False)},
        "paired_feature_comparisons": {
            "cox": {
                name: interval(v)
                for name, v in [("uno_c", 0.01), ("ibs", -0.01), ("mean_auc", 0.01)]
            }
        },
    }


def test_performance_exports_preserve_values_and_failed_ci(tmp_path):
    report = performance_fixture()
    report["models"]["cox__clinical"]["ibs"]["status"] = "conditional_approximate"
    original = deepcopy(report)
    paths = render_performance(report, tmp_path, "synthetic_performance", figure_config())
    assert report == original
    assert len(paths) == 19
    assert all(path.is_file() and path.stat().st_size > 100 for path in paths)
    table = (tmp_path / "tables" / "synthetic_performance_performance_table.tex").read_text()
    assert "0.710" in table and "0.720" in table
    assert "CI unavailable" in table
    assert "Integrated Brier score & 0.150 & 0.130 & 0.170" in table
    assert "intervals conditional on the defined bootstrap replicates" in table
    assert r"Cox PH$^{*}$" in table
    assert r"\toprule" in table and "S[table-format" in table
    paired = (tmp_path / "tables" / "synthetic_performance_paired_table.tex").read_text()
    assert "-0.010" in paired
    assert not plt.get_fignums()


def test_bayesian_exports_keep_gate_failure_and_prior_only_context(tmp_path):
    workflow = {
        "prior_predictive": {
            "time_histogram_edges_months": [1, 10, 100, 1000],
            "time_histogram_probability": [0.1, 0.5, 0.35],
            "extreme_fraction": 0.05,
        },
        "training_ppc": {
            "times_months": [12, 60, 120],
            "prior_survival_lower": [0.5, 0.3, 0.1],
            "prior_survival_upper": [0.99, 0.95, 0.8],
            "posterior_survival_lower": [0.9, 0.7, 0.5],
            "posterior_survival_upper": [0.95, 0.8, 0.6],
            "posterior_survival_mean": [0.93, 0.75, 0.55],
            "km_survival": [0.92, 0.73, 0.58],
            "km_lower": [0.88, 0.68, 0.52],
            "km_upper": [0.96, 0.79, 0.63],
        },
        "coefficients": {
            "age_at_diagnosis": {"mean_log_time_ratio": -0.2, "time_ratio_interval": [0.7, 0.95]},
            "receptor_unknown": {"mean_log_time_ratio": 0.0, "time_ratio_interval": [0.5, 2.0]},
        },
        "constant_training_features": ["receptor_unknown"],
    }
    diagnostics = {
        "max_rhat": 1.02,
        "min_ess_bulk": 900.0,
        "min_ess_tail": 750.0,
        "divergences": 0,
        "checks": {"rhat": False, "ess_bulk": True, "ess_tail": True, "divergences": True},
        "thresholds": {
            "max_rhat": 1.01,
            "min_ess_bulk": 400,
            "min_ess_tail": 400,
            "max_divergences": 0,
        },
    }
    original = deepcopy(workflow)
    paths = render_bayesian_workflow(
        workflow, diagnostics, tmp_path / "figures", "synthetic_bayesian", figure_config()
    )
    assert workflow == original
    assert len(paths) == 10
    assert all(path.is_file() and path.stat().st_size > 100 for path in paths)
    gate_table = (tmp_path / "tables" / "synthetic_bayesian_mcmc_diagnostics_table.tex").read_text()
    assert "1.0200" in gate_table and "FAIL" in gate_table
    coefficient_table = (
        tmp_path / "tables" / "synthetic_bayesian_coefficients_table.tex"
    ).read_text()
    assert f"{np.exp(-0.2):.3f}" in coefficient_table
    assert r"receptor\_unknown" in coefficient_table
    assert not plt.get_fignums()
