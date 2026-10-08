"""Rendering tests use aggregate synthetic summaries only, never patient outcomes."""

from copy import deepcopy

import pytest

from brca.viz.coverage import render_coverage


def synthetic_report():
    report = {}
    for model in ("cox", "weibull"):
        for alpha in (0.2, 0.1, 0.05):
            summary = {
                "alpha": alpha,
                "confidence_level": 0.95,
                "nominal_coverage": 1 - alpha,
                "horizon_months": 120,
                "n_cal": 200,
                "n_test": 80,
                "observable_coverage_lower": 0.72,
                "observable_coverage_upper": 0.96,
                "ipcw_coverage": 0.87,
                "ipcw_coverage_ci": [0.80, 0.94],
                "width_median_months": 85.3,
                "width_q25_months": 70.1,
                "width_q75_months": 99.2,
                "width_median_ci_months": [80.2, 90.4],
                "full_support_fraction": 0.10,
                "infinite_interval_fraction": 0.0,
            }
            subgroup = {**summary, "n_test": 22}
            report[f"{model}_{alpha}"] = {
                "endpoint": "disease_specific",
                "model": model,
                "feature_set": "clinical",
                "registered_method": "conservative",
                "marginal": summary,
                "subgroups": {
                    "pam50_subtype": {
                        "LumA": subgroup,
                        "rare_subtype": {
                            "n_test": 3,
                            "status": "suppressed_small_group",
                            "width_median_months": 987.6,
                            "ipcw_coverage": 0.123456,
                        },
                    }
                },
            }
    return report


def config():
    return {
        "figures": {
            "width_inches": 6.3,
            "height_inches": 4.0,
            "dpi": 90,
            "font_size": 9,
            "line_width": 1.6,
        },
        "conformal": {"alphas": [0.2, 0.1, 0.05], "horizons_months": [300, 120]},
    }


def test_renders_stage05_schema_without_recomputing_or_leaking_suppressed_values(tmp_path):
    report = synthetic_report()
    original = deepcopy(report)
    paths = render_coverage(report, tmp_path, config())
    assert report == original
    assert len(paths) == 7  # marginal/subgroups each PDF+PNG+table, one context caption
    assert all(path.is_file() and path.stat().st_size > 0 for path in paths)
    pdfs = [path for path in paths if path.suffix == ".pdf"]
    assert len(pdfs) == 2 and all(path.read_bytes().startswith(b"%PDF") for path in pdfs)
    assert len([path for path in paths if path.suffix == ".png"]) == 2
    tables = "\n".join(path.read_text() for path in paths if path.suffix == ".tex")
    assert "987.6" not in tables and "0.123456" not in tables
    assert "rare\\_subtype\\textsuperscript{s}" in tables
    assert "85.3" in tables and "[70.1, 99.2]" in tables
    assert "[0.80, 0.94]" in tables
    assert "\\toprule" in tables and "S[table-format" in tables
    caption = next(path for path in paths if path.suffix == ".txt").read_text()
    assert "not confidence intervals for latent" in caption
    assert "n_cal=200" in caption
    assert "Confidence level: 95%" in caption


def test_unavailable_ipcw_and_failed_interval_gate_remain_explicit(tmp_path):
    report = synthetic_report()
    for row in report.values():
        summary = row["marginal"]
        summary["ipcw_coverage"] = None
        summary["ipcw_coverage_ci"] = None
        summary["width_median_ci_months"] = [12.3, 45.6]
        summary["bootstrap_statistics"] = {
            "width_median_months": {"status": "bootstrap_gate_failed"}
        }
        row["subgroups"]["pam50_subtype"]["LumA"]["bootstrap_statistics"] = {
            "ipcw_coverage": {"status": "conditional_on_successful_resamples"}
        }
    paths = render_coverage(report, tmp_path, config())
    marginal_table = next(path for path in paths if "marginal" in path.name).read_text()
    assert "[12.3, 45.6]" not in marginal_table
    assert "85.3" in marginal_table
    subgroup_table = next(path for path in paths if "subgroups" in path.name).read_text()
    assert "[0.80, 0.94]\\textsuperscript{a}" in subgroup_table


def test_rejects_metrics_from_another_horizon_protocol(tmp_path):
    report = synthetic_report()
    next(iter(report.values()))["marginal"]["horizon_months"] = 240
    with pytest.raises(ValueError, match="alpha/horizon protocol"):
        render_coverage(report, tmp_path, config())
    assert not list(tmp_path.iterdir())
