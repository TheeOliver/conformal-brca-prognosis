"""Publication figures and tables rendered exclusively from saved aggregate metrics.

Pass the outputs root (or its figures directory). Returned paths include PDF,
PNG and booktabs/siunitx LaTeX so callers can record every artifact in a manifest.
No cohort, model, or metric-estimation code is imported here.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLOURS = ("#2a78d6", "#eb6834", "#1baf7a", "#6b4c9a", "#e87ba4")
REFERENCE = "#767676"
MODEL_NAMES = {
    "cox": "Cox PH",
    "coxph": "Cox PH",
    "cox_ph": "Cox PH",
    "rsf": "Random survival forest",
    "weibull": "Bayesian Weibull AFT",
    "bayesian": "Bayesian Weibull AFT",
    "bayesian_weibull": "Bayesian Weibull AFT",
    "bayesian_aft": "Bayesian Weibull AFT",
    "lognormal": "Bayesian lognormal AFT",
}
METRIC_NAMES = {
    "uno_c": "Uno C-index",
    "ibs": "Integrated Brier score",
    "mean_auc": "Mean time-dependent AUC",
    "brier": "Brier score",
    "auc": "Time-dependent AUC",
}


def _directories(output_dir: Path | str) -> tuple[Path, Path]:
    root = Path(output_dir)
    figures = root if root.name == "figures" else root / "figures"
    tables = figures.parent / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    return figures, tables


def _style(config) -> None:
    plt.rcParams.update(
        {
            "font.size": max(8, config.figures.font_size),
            "axes.spines.top": False,
            "axes.spines.right": False,
            "lines.solid_capstyle": "round",
            "pdf.fonttype": 42,
        }
    )


def _save(fig, directory: Path, name: str, config) -> list[Path]:
    paths = [directory / f"{name}.pdf", directory / f"{name}.png"]
    # Keep the declared physical width. Cropping to a tight bounding box can
    # silently make long coefficient labels exceed the thesis text width.
    fig.set_layout_engine("constrained")
    fig.savefig(paths[0], metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(paths[1], dpi=config.figures.dpi)
    plt.close(fig)
    return paths


def _escape(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in str(value))


def _number(value, precision=3) -> str:
    return (
        f"{float(value):.{precision}f}"
        if value is not None and np.isfinite(value)
        else r"\multicolumn{1}{c}{--}"
    )


def _write_table(
    path: Path, columns: str, headers: list[str], rows: list[list[str]], note: str
) -> Path:
    # S columns require siunitx, avoiding fake alignment with spaces or differing
    # decimal precision. The table can be input directly into the thesis.
    lines = [
        "% Requires \\usepackage{booktabs,siunitx}",
        f"% {note}",
        r"\begingroup\footnotesize\setlength{\tabcolsep}{3pt}",
        f"\\begin{{tabular}}{{{columns}}}",
        r"\toprule",
        " & ".join(headers) + r" \\",
        r"\midrule",
    ]
    lines.extend(" & ".join(row) + r" \\" for row in rows)
    if any(r"$^{*}$" in cell for row in rows for cell in row):
        note += " * marks intervals conditional on the defined bootstrap replicates."
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\par\smallskip{\footnotesize " + _escape(note) + "}",
            r"\endgroup",
            "",
        ]
    )
    path.write_text("\n".join(lines))
    return path


def _labels(identifier: str) -> tuple[str, str]:
    model, _, arm = identifier.partition("__")
    name = MODEL_NAMES.get(model, model.replace("_", " "))
    feature = {"clinical": "Clinical", "clinical_molecular": "Clinical + PAM50"}.get(
        arm, arm.replace("_", " ")
    )
    return name, feature


def _has_ci(summary: dict) -> bool:
    values = [summary.get("ci_lower"), summary.get("ci_upper")]
    return summary.get("status", "ok") in {"ok", "conditional_approximate"} and all(
        value is not None and np.isfinite(value) for value in values
    )


def _table_name(name: str, summary: dict) -> str:
    suffix = r"$^{*}$" if summary.get("status") == "conditional_approximate" else ""
    return _escape(name) + suffix


def _finite_estimate(summary: dict) -> bool:
    value = summary.get("estimate")
    return value is not None and np.isfinite(value)


def _table_interval(summary: dict, precision=3) -> list[str]:
    if _has_ci(summary):
        return [_number(summary["ci_lower"], precision), _number(summary["ci_upper"], precision)]
    return [r"\multicolumn{2}{c}{CI unavailable}"]


def _scalar_plot(
    entries: dict, metric: str, config, confidence: float, *, difference: bool = False
):
    direction = metric != "ibs"
    keys = sorted(
        entries,
        key=lambda key: (
            entries[key][metric]["estimate"] if _finite_estimate(entries[key][metric]) else -np.inf
        ),
        reverse=direction,
    )
    fig, ax = plt.subplots(
        figsize=(config.figures.width_inches, max(config.figures.height_inches, len(keys) * 0.42))
    )
    labels = []
    for position, key in enumerate(keys):
        summary = entries[key][metric]
        model, arm = _labels(key)
        label = model if difference else f"{model}\n{arm}"
        molecular = key.endswith("__clinical_molecular")
        colour, marker = (COLOURS[1], "D") if molecular else (COLOURS[0], "o")
        if _has_ci(summary):
            ax.hlines(
                position,
                summary["ci_lower"],
                summary["ci_upper"],
                color=colour,
                linewidth=max(2, config.figures.line_width),
            )
        else:
            label += " (CI unavailable)"
        if summary.get("status") == "conditional_approximate":
            label += " (conditional CI)"
        if _finite_estimate(summary):
            ax.plot(
                summary["estimate"],
                position,
                marker=marker,
                markersize=8,
                markerfacecolor=colour if _has_ci(summary) else "white",
                markeredgecolor=colour,
                linestyle="none",
            )
        else:
            label += " (estimate unavailable)"
        labels.append(label)
    ax.set_yticks(range(len(keys)), labels)
    ax.invert_yaxis()
    if difference:
        ax.axvline(0, color=REFERENCE, linewidth=1)
        suffix = "positive favours PAM50" if direction else "negative favours PAM50"
        ax.set_xlabel(
            f"Difference in {METRIC_NAMES[metric]}\nClinical + PAM50 minus clinical ({suffix})"
        )
    else:
        better = "higher is better" if direction else "lower is better"
        ax.set_xlabel(f"{METRIC_NAMES[metric]} ({better}; unitless)")
    ax.set_xlabel(ax.get_xlabel() + f"\n{confidence:.0%} bootstrap intervals")
    return fig


def render_performance(report: dict, output_dir: Path | str, stem: str, config) -> list[Path]:
    """Export fixed-model performance and paired feature differences from evaluate_models."""
    figures, tables = _directories(output_dir)
    _style(config)
    models = report["models"]
    if not models:
        raise ValueError("Performance report contains no model summaries")
    paths = []
    confidence = report["bootstrap"]["confidence_level"]
    scalar_rows = []
    curve_rows = []
    for metric in ("uno_c", "ibs", "mean_auc"):
        paths.extend(
            _save(
                _scalar_plot(models, metric, config, confidence),
                figures,
                f"{stem}_{metric}",
                config,
            )
        )
        for key, summaries in sorted(models.items()):
            summary = summaries[metric]
            model, arm = _labels(key)
            scalar_rows.append(
                [
                    _table_name(model, summary),
                    _escape(arm),
                    _escape(METRIC_NAMES[metric]),
                    _number(summary.get("estimate")),
                    *_table_interval(summary),
                ]
            )
    paths.append(
        _write_table(
            tables / f"{stem}_performance_table.tex",
            "lll S[table-format=1.3] S[table-format=1.3] S[table-format=1.3]",
            ["Model", "Features", "Metric", "{Estimate}", "{CI lower}", "{CI upper}"],
            scalar_rows,
            f"{confidence:.0%} paired test-patient bootstrap intervals; fixed fits and "
            "training censoring weights. CI unavailable means the bootstrap gate failed.",
        )
    )
    base_models = sorted({key.partition("__")[0] for key in models})
    colours = {name: COLOURS[j % len(COLOURS)] for j, name in enumerate(base_models)}
    markers = {name: ("o", "s", "^")[j % 3] for j, name in enumerate(base_models)}
    for metric in ("brier", "auc"):
        fig, ax = plt.subplots(figsize=(config.figures.width_inches, config.figures.height_inches))
        for key, summaries in sorted(models.items()):
            points = summaries[metric]
            time = np.array([point["time_months"] for point in points])
            value = np.array([point.get("estimate") for point in points], dtype=float)
            model, arm = _labels(key)
            colour = colours[key.partition("__")[0]]
            molecular = key.endswith("__clinical_molecular")
            label = f"{model}: {arm}"
            if not all(_has_ci(point) for point in points):
                label += " (some CI unavailable)"
            if any(point.get("status") == "conditional_approximate" for point in points):
                label += " (conditional CI)"
            ax.plot(
                time,
                value,
                color=colour,
                linestyle="--" if molecular else "-",
                marker=markers[key.partition("__")[0]],
                markersize=8,
                linewidth=max(2, config.figures.line_width),
                label=label,
            )
            lower = np.array([point["ci_lower"] if _has_ci(point) else np.nan for point in points])
            upper = np.array([point["ci_upper"] if _has_ci(point) else np.nan for point in points])
            ax.fill_between(time, lower, upper, color=colour, alpha=0.10)
            for point in points:
                curve_rows.append(
                    [
                        _table_name(model, point),
                        _escape(arm),
                        _escape(METRIC_NAMES[metric]),
                        _number(point["time_months"], 1),
                        _number(point.get("estimate")),
                        *_table_interval(point),
                    ]
                )
        ax.set(xlabel="Time (months)", ylabel=f"{METRIC_NAMES[metric]} (unitless)")
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.20),
            fontsize=max(8, config.figures.font_size),
            frameon=False,
        )
        paths.extend(_save(fig, figures, f"{stem}_{metric}", config))
    paths.append(
        _write_table(
            tables / f"{stem}_time_metrics_table.tex",
            "lll S[table-format=3.1] S[table-format=1.3] S[table-format=1.3] S[table-format=1.3]",
            ["Model", "Features", "Metric", "{Months}", "{Estimate}", "{CI lower}", "{CI upper}"],
            curve_rows,
            f"Pointwise {confidence:.0%} bootstrap intervals on the fixed "
            "training-derived time grid.",
        )
    )
    comparisons = report.get("paired_feature_comparisons", {})
    if comparisons:
        paired_rows = []
        for metric in ("uno_c", "ibs", "mean_auc"):
            paths.extend(
                _save(
                    _scalar_plot(comparisons, metric, config, confidence, difference=True),
                    figures,
                    f"{stem}_paired_{metric}",
                    config,
                )
            )
            for key, summary in sorted(comparisons.items()):
                paired_rows.append(
                    [
                        _table_name(_labels(key)[0], summary[metric]),
                        _escape(METRIC_NAMES[metric]),
                        _number(summary[metric].get("estimate"), 5),
                        *_table_interval(summary[metric], 5),
                    ]
                )
        paths.append(
            _write_table(
                tables / f"{stem}_paired_table.tex",
                "ll S[table-format=-1.5] S[table-format=-1.5] S[table-format=-1.5]",
                ["Model", "Metric", "{Difference}", "{CI lower}", "{CI upper}"],
                paired_rows,
                f"Clinical + PAM50 minus clinical; {confidence:.0%} paired bootstrap intervals. "
                "Positive favours molecular for C/AUC; negative for IBS. "
                "No multiplicity adjustment.",
            )
        )
    return paths


def render_bayesian_workflow(
    workflow: dict, diagnostics: dict, output_dir: Path | str, stem: str, config
) -> list[Path]:
    """Render prior checks, training PPC, coefficients and recorded MCMC diagnostics."""
    figures, tables = _directories(output_dir)
    _style(config)
    paths = []
    prior = workflow["prior_predictive"]
    fig, ax = plt.subplots(figsize=(config.figures.width_inches, config.figures.height_inches))
    edges = np.asarray(prior["time_histogram_edges_months"])
    ax.bar(
        edges[:-1],
        prior["time_histogram_probability"],
        width=np.diff(edges),
        align="edge",
        color=COLOURS[0],
        edgecolor="white",
    )
    ax.set(
        xscale="log",
        xlabel="Prior event time (months, log scale)",
        ylabel="Prior probability per log-time bin",
    )
    ax.text(
        0.02,
        0.98,
        f"Mass outside displayed bounds: {prior['extreme_fraction']:.1%}",
        va="top",
        transform=ax.transAxes,
    )
    paths.extend(_save(fig, figures, f"{stem}_prior_times", config))
    ppc = workflow.get("training_ppc")
    probability = config.bayesian.predictive_interval_probability
    if ppc:
        fig, ax = plt.subplots(figsize=(config.figures.width_inches, config.figures.height_inches))
        time = np.array(ppc["times_months"])
        ax.fill_between(
            time,
            ppc["prior_survival_lower"],
            ppc["prior_survival_upper"],
            color=COLOURS[0],
            alpha=0.15,
            hatch="//",
            label=f"Prior {probability:.0%} pointwise band",
        )
        ax.fill_between(
            time,
            ppc["posterior_survival_lower"],
            ppc["posterior_survival_upper"],
            color=COLOURS[1],
            alpha=0.20,
            label=f"Posterior {probability:.0%} pointwise band",
        )
        ax.plot(
            time,
            ppc["posterior_survival_mean"],
            color=COLOURS[1],
            linestyle="-",
            linewidth=max(2, config.figures.line_width),
            label="Posterior mean survival",
        )
        # Only grid evaluations are present in the aggregate; do not invent a
        # full Kaplan-Meier step curve with jumps at the grid rather than events.
        ax.vlines(
            time,
            ppc["km_lower"],
            ppc["km_upper"],
            color=COLOURS[2],
            linewidth=max(2, config.figures.line_width),
        )
        ax.plot(
            time,
            ppc["km_survival"],
            marker="s",
            markersize=8,
            linestyle="none",
            color=COLOURS[2],
            label=f"Training KM at grid with {probability:.0%} CI",
        )
        ax.set(xlabel="Time (months)", ylabel="Survival probability", ylim=(0, 1))
        ax.legend(frameon=False, fontsize=max(8, config.figures.font_size))
        paths.extend(_save(fig, figures, f"{stem}_training_ppc", config))
    coefficients = workflow.get("coefficients", {})
    if coefficients:
        names = list(coefficients)
        fig, ax = plt.subplots(
            figsize=(
                config.figures.width_inches,
                max(config.figures.height_inches, len(names) * 0.26),
            )
        )
        rows = []
        labels = []
        for j, name in enumerate(names):
            item = coefficients[name]
            lower, upper = item["time_ratio_interval"]
            point = np.exp(item["mean_log_time_ratio"])
            ax.hlines(
                j, lower, upper, color=COLOURS[0], linewidth=max(2, config.figures.line_width)
            )
            ax.plot(point, j, "o", color=COLOURS[0], markersize=8)
            labels.append(
                name.replace("_", " ")
                + (
                    " [constant in train]"
                    if name in workflow.get("constant_training_features", [])
                    else ""
                )
            )
            rows.append([_escape(name), _number(point), _number(lower), _number(upper)])
        ax.set_yticks(range(len(names)), labels)
        ax.invert_yaxis()
        ax.axvline(1, color=REFERENCE, linewidth=1)
        ax.set(
            xscale="log",
            xlabel=f"Conditional time ratio (log scale)\n{probability:.0%} credible interval",
        )
        paths.extend(_save(fig, figures, f"{stem}_coefficients", config))
        paths.append(
            _write_table(
                tables / f"{stem}_coefficients_table.tex",
                "l S[table-format=3.3] S[table-format=3.3] S[table-format=3.3]",
                ["Predictor", "{Geometric mean}", "{CrI lower}", "{CrI upper}"],
                rows,
                "Point is exp(posterior mean log time ratio); associations conditional on other "
                "predictors. Continuous predictors are per training standard deviation. "
                "Intervals use all retained posterior draws.",
            )
        )
    threshold = diagnostics["thresholds"]
    diagnostic_rows = [
        [
            "Maximum R-hat",
            f"{diagnostics['max_rhat']:.4f}",
            f"<= {threshold['max_rhat']}",
            "pass" if diagnostics["checks"]["rhat"] else "FAIL",
        ],
        [
            "Minimum bulk ESS",
            f"{diagnostics['min_ess_bulk']:.1f}",
            f">= {threshold['min_ess_bulk']}",
            "pass" if diagnostics["checks"]["ess_bulk"] else "FAIL",
        ],
        [
            "Minimum tail ESS",
            f"{diagnostics['min_ess_tail']:.1f}",
            f">= {threshold['min_ess_tail']}",
            "pass" if diagnostics["checks"]["ess_tail"] else "FAIL",
        ],
        [
            "Divergences",
            str(diagnostics["divergences"]),
            f"<= {threshold['max_divergences']}",
            "pass" if diagnostics["checks"]["divergences"] else "FAIL",
        ],
    ]
    fig, ax = plt.subplots(figsize=(config.figures.width_inches, config.figures.height_inches))
    ax.axis("off")
    table = ax.table(
        cellText=diagnostic_rows,
        colLabels=["Diagnostic", "Observed", "Required", "Gate"],
        cellLoc="left",
        loc="center",
        colWidths=[0.38, 0.20, 0.26, 0.16],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(max(8, config.figures.font_size))
    table.scale(1, 1.8)
    for cell in table.get_celld().values():
        cell.set_linewidth(0)
    paths.extend(_save(fig, figures, f"{stem}_mcmc_diagnostics", config))
    rows = [
        [
            _escape(row[0]),
            _number(float(row[1]), 4 if "R-hat" in row[0] else 1),
            _escape(row[2]),
            row[3],
        ]
        for row in diagnostic_rows
    ]
    paths.append(
        _write_table(
            tables / f"{stem}_mcmc_diagnostics_table.tex",
            "l S[table-format=5.4] ll",
            ["Diagnostic", "{Observed}", "Required", "Gate"],
            rows,
            "Gates are copied from the saved diagnostics, not recomputed by the renderer.",
        )
    )
    return paths
