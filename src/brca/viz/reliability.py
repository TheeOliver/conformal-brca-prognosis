"""Export saved uncertainty-reliability comparisons; never recompute statistics."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter

METHODS = {
    "conservative_score_envelope": ("Primary two-sided", "#2a78d6", "o", -0.27),
    "conservative_lower_bound": ("One-sided lower", "#eb6834", "s", -0.09),
    "qin_bootstrap_cox": ("Original Qin Cox", "#6b4c9a", "D", 0.09),
    "qin_bootstrap_ridge_cox": ("Ridge Qin Cox", "#1baf7a", "^", 0.27),
}
MODELS = {"cox": "Cox", "cox_ridge": "Cox ridge", "weibull": "Weibull", "rsf": "RSF"}
FEATURES = {"clinical": "Clinical", "clinical_molecular": "+ PAM50"}
ENDPOINTS = {"disease_specific": "Breast-cancer-specific survival", "overall": "Overall survival"}


def _finite(value):
    return value is not None and isinstance(value, (int, float)) and np.isfinite(value)


def _fmt(value, *, scale=1.0, digits=1, signed=False):
    if not _finite(value):
        return "unavailable"
    return format(value * scale, f"{'+' if signed else ''}.{digits}f")


def _usable_ci(statistic):
    status = statistic.get("status", "ok")
    if status == "bootstrap_gate_failed" or status.startswith("unavailable"):
        return None
    ci = statistic.get("ci")
    return ci if ci is not None and len(ci) == 2 and all(_finite(x) for x in ci) else None


def _estimate_ci(estimate, statistic, *, scale=1.0, digits=1, signed=False):
    if not _finite(estimate):
        return "unavailable"
    text = _fmt(estimate, scale=scale, digits=digits, signed=signed)
    ci = _usable_ci(statistic)
    if ci is not None:
        text += " [" + ", ".join(_fmt(x, scale=scale, digits=digits) for x in ci) + "]"
    else:
        text += " [CI unavailable]"
    if statistic.get("status") == "conditional_on_successful_resamples":
        text += "†"
    return text


def _summary_stat(summary, name, ci_name):
    return {"ci": summary.get(ci_name), **summary.get("bootstrap_statistics", {}).get(name, {})}


def _difference(comparison, name, *, scale=1.0):
    item = comparison["differences"][name]
    return _estimate_ci(item.get("estimate"), item, scale=scale, signed=True)


def _method(row):
    return row["marginal"]["method"]


def _context(row):
    summary = row["marginal"]
    return row["endpoint"], summary["horizon_months"], summary["alpha"]


def _sort_row(row):
    model = "cox" if row["model"] == "cox_ridge" else row["model"]
    return (
        list(MODELS).index(model),
        list(FEATURES).index(row["feature_set"]),
        list(METHODS).index(_method(row)),
    )


def _all_rows(inputs):
    rows = {
        key: row
        for key, row in inputs["primary_coverage"].items()
        if _method(row) == "conservative_score_envelope"
    }
    rows.update(inputs["original_qin_coverage"])
    rows.update(inputs["coverage"])
    for row in rows.values():
        if _method(row) not in METHODS or row["model"] not in MODELS:
            raise ValueError("Reliability export needs a registered method/model style")
    return rows


def _coverage_table(rows):
    lines = [
        "| Model | Features | Method | n_cal / B | Observable coverage bounds (%) | "
        "IPCW coverage (%) [CI] | Median width (months) [CI] | Width IQR (months) |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted(rows, key=_sort_row):
        s = row["marginal"]
        bounds = " to ".join(
            _fmt(s.get(key), scale=100)
            for key in ("observable_coverage_lower", "observable_coverage_upper")
        )
        ipcw = _estimate_ci(
            s.get("ipcw_coverage"), _summary_stat(s, "ipcw_coverage", "ipcw_coverage_ci"), scale=100
        )
        width = _estimate_ci(
            s.get("width_median_months"),
            _summary_stat(s, "width_median_months", "width_median_ci_months"),
        )
        iqr = " to ".join(_fmt(s.get(key)) for key in ("width_q25_months", "width_q75_months"))
        lines.append(
            f"| {MODELS[row['model']]} | {FEATURES[row['feature_set']]} | "
            f"{METHODS[_method(row)][0]} | {s['n_cal']} | {bounds} | {ipcw} | {width} | {iqr} |"
        )
    return lines


def _paired_table(pairs, coverage, context, *, features=False):
    lines = [
        "| Model | Features / candidate | Δ observable lower (pp) [CI] | "
        "Δ observable upper (pp) [CI] | Δ IPCW coverage (pp) [CI] | "
        "Δ median width (months) [CI] |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    keys = sorted(
        (key for key in pairs if _context(coverage[key]) == context),
        key=lambda key: _sort_row(coverage[key]),
    )
    for key in keys:
        row, pair = coverage[key], pairs[key]
        label = METHODS[_method(row)][0]
        if not features:
            label = f"{FEATURES[row['feature_set']]} / {label}"
        lines.append(
            f"| {MODELS[row['model']]} | {label} | "
            f"{_difference(pair, 'observable_coverage_lower', scale=100)} | "
            f"{_difference(pair, 'observable_coverage_upper', scale=100)} | "
            f"{_difference(pair, 'ipcw_coverage', scale=100)} | "
            f"{_difference(pair, 'width_median_months')} |"
        )
    return lines


def _diagnostic_tables(inputs):
    preparation = inputs["preparation"]
    lines = [
        "## Training diagnostics",
        "",
        "Ridge and original Qin use the same training event pool. The unpenalized pilot "
        "is a shorter diagnostic run; its denominator differs from the complete original run.",
        "",
        "| Endpoint | Features | Original Qin refits | Pilot refits | Ridge refits | "
        "Ridge resamples with constant columns | Event-pool effective size |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for _key, item in sorted(preparation["ridge"].items()):
        original = next(
            row
            for row in inputs["original_qin_coverage"].values()
            if row["endpoint"] == item["endpoint"] and row["feature_set"] == item["feature_set"]
        )
        old, new = original["bootstrap_diagnostics"], item["diagnostics"]
        pilot = item["unpenalized_pilot"]
        lines.append(
            f"| {ENDPOINTS[item['endpoint']]} | {FEATURES[item['feature_set']]} | "
            f"{old['bootstrap_successful']}/{old['bootstrap_requested']} | "
            f"{pilot['counts'].get('successful', 0)}/{pilot['requested']} | "
            f"{new['bootstrap_successful']}/{new['bootstrap_requested']} | "
            f"{new['bootstrap_design']['replicates_with_constant_columns']} | "
            f"{_fmt(new['event_pool_effective_sample_size'])} |"
        )
    lines += [
        "",
        "| Endpoint | Horizon (months) | Training G(tau-) | IPCW support | "
        "Training patients observed through horizon | Restricted known-target effective size |",
        "| --- | ---: | ---: | --- | ---: | ---: |",
    ]
    for endpoint, support in sorted(preparation["training_support"].items()):
        for item in sorted(support["horizons"].values(), key=lambda item: item["horizon_months"]):
            lines.append(
                f"| {ENDPOINTS[endpoint]} | {item['horizon_months']:g} | "
                f"{_fmt(item['censoring_survival_at_horizon'], digits=4)} | "
                f"{item['ipcw_support_status']} | {item['n_observed_through_horizon']} | "
                f"{_fmt(item['restricted_known_target_effective_sample_size'])} |"
            )
    return lines


def _performance_tables(inputs):
    lines = [
        "## Original and ridge Cox probability/ranking metrics",
        "",
        "These are separately stored estimates and confidence intervals. No paired "
        "original-versus-ridge performance difference was computed; overlap of marginal "
        "intervals is not a test of a difference. Higher Uno C is better; lower IBS is better.",
        "",
        "| Endpoint | Features | Cox fit | Uno C [CI] | IBS [CI] | Integration grid (months) |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for endpoint, ridge in sorted(inputs["ridge_performance"].items()):
        for feature in FEATURES:
            for label, report, key in (
                ("Original", inputs["primary_performance"][endpoint], f"cox__{feature}"),
                ("Ridge", ridge, f"cox_ridge__{feature}"),
            ):
                row = report["models"][key]
                values = []
                for metric in ("uno_c", "ibs"):
                    statistic = row[metric]
                    ci = {**statistic, "ci": [statistic.get("ci_lower"), statistic.get("ci_upper")]}
                    values.append(_estimate_ci(statistic.get("estimate"), ci, digits=3))
                grid = report["time_grid_months"]
                lines.append(
                    f"| {ENDPOINTS[endpoint]} | {FEATURES[feature]} | {label} | "
                    f"{values[0]} | {values[1]} | {grid[0]:.3f} to {grid[-1]:.3f} |"
                )
    return lines


def _report(inputs, rows, figure_paths):
    confidence = next(iter(inputs["coverage"].values()))["marginal"]["confidence_level"]
    n_test = sorted({row["marginal"]["n_test"] for row in rows.values()})
    lines = [
        "# Exploratory survival-uncertainty reliability comparison",
        "",
        "All numbers below are read from saved aggregate evaluation outputs. The original "
        "test set has already been used: these comparisons are post hoc and cannot select "
        "a newly validated primary method. Frozen primary results remain the reference.",
        "",
        "The one-sided method reports **[L, τ] for min(T, τ)**, equivalently a lower "
        "survival bound T ≥ L. The upper endpoint τ is a target boundary, not a predicted "
        "time by which survival ends. Primary and Qin intervals answer a two-sided "
        "question. Ridge addresses numerical stability; it does not repair censoring "
        "assumptions or establish clinical utility.",
        "",
        f"Test sample size(s): {', '.join(map(str, n_test))}. Confidence intervals are "
        f"{confidence:.0%} paired test-patient bootstrap intervals conditional on fitted "
        "models, calibration, and the training censoring fit. Observable coverage bounds "
        "are identification bounds, **not confidence intervals**. IPCW is assumption-dependent "
        "and shown only where supported. † marks intervals conditional on successful "
        "bootstrap resamples. Unavailable values are never replaced by zero.",
        "",
        "For primary/lower methods, n_cal counts held-out calibration patients. For Qin "
        "methods, B counts successful bootstrap pivots; it is not a held-out sample size. "
        "All widths and width differences are months. pp denotes percentage points.",
        "",
        *_diagnostic_tables(inputs),
        "",
        *_performance_tables(inputs),
        "",
        "## Coverage, width, and paired comparisons at every registered setting",
        "",
        "Method differences are **one-sided lower minus primary two-sided**, or "
        "**ridge Qin minus original Qin**, for the same model/feature arm. Positive "
        "width differences mean wider candidate intervals. Feature differences are "
        "**clinical + PAM50 minus clinical** within the same candidate method. "
        "These intervals quantify test resampling uncertainty, not full retraining uncertainty.",
    ]
    contexts = sorted({_context(row) for row in rows.values()})
    for context in contexts:
        endpoint, horizon, alpha = context
        selected = [row for row in rows.values() if _context(row) == context]
        lines += [
            "",
            f"### {ENDPOINTS[endpoint]} · {horizon:g} months · {1 - alpha:.0%} nominal",
            "",
            *_coverage_table(selected),
            "",
            "Paired method differences:",
            "",
            *_paired_table(inputs["paired_method_comparisons"], inputs["coverage"], context),
            "",
            "Paired feature differences:",
            "",
            *_paired_table(
                inputs["paired_feature_comparisons"], inputs["coverage"], context, features=True
            ),
        ]
    lines += ["", "## Figures and provenance", ""]
    for path in figure_paths:
        if path.suffix == ".pdf":
            preview = path.with_suffix(".png")
            lines.append(
                f"- [{path.name}](../figures/{path.name}) "
                f"([PNG preview](../figures/{preview.name}))"
            )
    lines += [
        "",
        "Suggested caption: Exploratory comparison on the reused test cohort. The top "
        "panel shows observable coverage bounds (thick ranges), IPCW estimates (filled "
        "symbols) with bootstrap confidence intervals (thin whiskers), and nominal coverage "
        "(grey line). Hollow symbols mark observable bounds where IPCW is unavailable. "
        "The lower panel shows median restricted interval width, IQR (thick range), and "
        "bootstrap median confidence interval (thin whiskers). One-sided lower bounds "
        "retain the horizon as a fixed upper boundary and answer a different question "
        "from the two-sided methods. Clinical and +PAM50 use identical patients.",
        "",
        "Every source file SHA-256 and the supporting evaluation manifest are recorded "
        "in `export_metadata.json`; the export also has its own StageRun manifest. "
        "See `coverage.json` for observable-bound confidence intervals and all subgroup "
        "summaries. Small-group suppression follows the original analysis protocol.",
        "",
    ]
    return "\n".join(lines)


def _plot_summary(axes, x, row):
    s = row["marginal"]
    _, color, marker, _ = METHODS[_method(row)]
    lo, hi = s["observable_coverage_lower"], s["observable_coverage_upper"]
    axes[0].vlines(x, lo, hi, color=color, linewidth=5, alpha=0.3)
    ipcw = s.get("ipcw_coverage")
    if _finite(ipcw):
        ci = _usable_ci(_summary_stat(s, "ipcw_coverage", "ipcw_coverage_ci"))
        if ci is not None:
            axes[0].vlines(x, *ci, color=color, linewidth=1.2)
        axes[0].plot(
            x,
            ipcw,
            marker=marker,
            color=color,
            markersize=8,
            markeredgecolor="white",
            markeredgewidth=0.8,
            linestyle="none",
        )
    else:
        axes[0].plot(
            [x, x],
            [lo, hi],
            marker=marker,
            color=color,
            markersize=8,
            markerfacecolor="none",
            linestyle="none",
        )
    axes[1].vlines(
        x, s["width_q25_months"], s["width_q75_months"], color=color, linewidth=5, alpha=0.3
    )
    ci = _usable_ci(_summary_stat(s, "width_median_months", "width_median_ci_months"))
    if ci is not None:
        axes[1].vlines(x, *ci, color=color, linewidth=1.2)
    axes[1].plot(
        x,
        s["width_median_months"],
        marker=marker,
        color=color,
        markersize=8,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linestyle="none",
    )


def _figure(rows, context, output_dir, settings):
    endpoint, horizon, alpha = context
    positions = [(model, feature) for model in ("cox", "weibull", "rsf") for feature in FEATURES]
    fig, axes = plt.subplots(
        2, 1, sharex=True, figsize=(settings.width_inches, settings.height_inches * 1.8)
    )
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.13, top=0.77, hspace=0.18)
    coverage_values = []
    for row in rows:
        model = "cox" if row["model"] == "cox_ridge" else row["model"]
        x = positions.index((model, row["feature_set"])) + METHODS[_method(row)][3]
        _plot_summary(axes, x, row)
        s = row["marginal"]
        coverage_values += [s["observable_coverage_lower"], s["observable_coverage_upper"]]
        ci = _usable_ci(_summary_stat(s, "ipcw_coverage", "ipcw_coverage_ci"))
        if ci is not None:
            coverage_values += ci
    axes[0].axhline(1 - alpha, color="#767676", linewidth=1)
    axes[0].set(ylabel="Coverage (%)", ylim=(max(0, min(coverage_values) - 0.05), 1.02))
    axes[0].yaxis.set_major_formatter(PercentFormatter(1))
    axes[1].set(ylabel="Interval width (months)", ylim=(0, horizon * 1.04))
    axes[1].set_xticks(
        range(len(positions)),
        [f"{MODELS[model]}\n{FEATURES[feature]}" for model, feature in positions],
    )
    axes[1].set_xlim(-0.6, len(positions) - 0.4)
    for axis in axes:
        axis.grid(axis="y", color="#dddddd", linewidth=0.5)
        axis.set_axisbelow(True)
    fig.legend(
        handles=[
            Line2D([], [], color=color, marker=marker, linestyle="none", markersize=8, label=label)
            for label, color, marker, _ in METHODS.values()
        ],
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.55, 0.995),
    )
    fig.text(
        0.14,
        0.835,
        "Coverage: thick bounds; filled IPCW / thin CI.\n"
        "Hollow bounds: IPCW unavailable. Width: median / IQR / CI.",
        fontsize=8,
    )
    fig.text(
        0.14,
        0.785,
        "One-sided lower uses a fixed horizon upper edge; test cohort reused.",
        fontsize=8,
    )
    n_test = sorted({row["marginal"]["n_test"] for row in rows})
    n_cal = sorted(
        {
            row["marginal"]["n_cal"]
            for row in rows
            if _method(row) in ("conservative_score_envelope", "conservative_lower_bound")
        }
    )
    n_pivots = sorted(
        {
            row["marginal"]["n_cal"]
            for row in rows
            if _method(row) in ("qin_bootstrap_cox", "qin_bootstrap_ridge_cox")
        }
    )
    counts = "; ".join(
        (
            "n_test=" + ",".join(map(str, n_test)),
            "primary/lower n_cal=" + ",".join(map(str, n_cal)),
            "Qin B=" + ",".join(map(str, n_pivots)),
        )
    )
    confidence = rows[0]["marginal"]["confidence_level"]
    fig.text(
        0.14, 0.035, f"{counts}\nWhiskers: {confidence:.0%} fixed-fit bootstrap CI.", fontsize=8
    )
    stem = f"reliability_coverage_width__{endpoint}__horizon{horizon:g}__nominal{1 - alpha:.2f}"
    paths = [output_dir / f"{stem}.pdf", output_dir / f"{stem}.png"]
    fig.savefig(paths[0], metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(paths[1], dpi=settings.dpi)
    plt.close(fig)
    return paths


def render_reliability(inputs: dict, output_dir: Path | str, config) -> list[Path]:
    """Render all supplied aggregate settings and four registered-primary-alpha figures."""
    output_dir = Path(output_dir)
    figure_dir, report_dir = output_dir / "figures", output_dir / "uncertainty_reliability"
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    rows = _all_rows(inputs)
    alpha = inputs["plan"]["settings"]["report_alpha"]
    contexts = sorted({_context(row) for row in rows.values() if row["marginal"]["alpha"] == alpha})
    paths = []
    with plt.rc_context(
        {
            "font.size": max(8, config.figures.font_size),
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    ):
        for context in contexts:
            selected = [row for row in rows.values() if _context(row) == context]
            paths += _figure(selected, context, figure_dir, config.figures)
    report = report_dir / "comparison_report.md"
    report.write_text(_report(inputs, rows, paths), encoding="utf-8")
    return [*paths, report]
