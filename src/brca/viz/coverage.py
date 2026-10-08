"""Publication exports of stored stage-05 coverage summaries, without recomputing metrics."""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

MODEL_STYLE = {
    "cox": ("Cox PH", "#2a78d6", "o", "-"),
    "weibull": ("Bayesian Weibull", "#eb6834", "s", "--"),
    "rsf": ("Random survival forest", "#1baf7a", "^", "-."),
    "lognormal": ("Bayesian log-normal", "#6b4c9a", "D", ":"),
}
REFERENCE = "#767676"


def _get(obj: Any, key: str) -> Any:
    return obj[key] if isinstance(obj, Mapping) else getattr(obj, key)


def _slug(value: Any) -> str:
    return re.sub(r"[^A-Za-z0-9-]+", "_", str(value)).strip("_").lower()


def _tex(value: Any) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "_": r"\_",
        "%": r"\%",
        "&": r"\&",
        "#": r"\#",
        "$": r"\$",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(char, char) for char in str(value))


def _style(model: str) -> tuple[str, str, str, str]:
    if model not in MODEL_STYLE:
        raise ValueError(f"Define a stable publication style for model {model!r}.")
    return MODEL_STYLE[model]


def _number(value: Any) -> bool:
    return value is not None and isinstance(value, (int, float)) and np.isfinite(value)


def _visible(summary: dict) -> bool:
    return summary.get("status") != "suppressed_small_group" and summary.get("n_test", 0) > 0


def _ci(summary: dict, statistic: str, key: str) -> list | None:
    status = summary.get("bootstrap_statistics", {}).get(statistic, {}).get("status", "ok")
    if status == "bootstrap_gate_failed" or status.startswith("unavailable"):
        return None
    value = summary.get(key)
    return (
        value if value is not None and len(value) == 2 and all(_number(x) for x in value) else None
    )


def _partial_ci(summary: dict, statistic: str) -> bool:
    return summary.get("bootstrap_statistics", {}).get(statistic, {}).get("status") == (
        "conditional_on_successful_resamples"
    )


def _context_name(context: tuple) -> str:
    endpoint, features, method, horizon = context
    return "__".join(map(_slug, (endpoint, features, method, f"horizon{horizon:g}")))


def _save(fig, directory: Path, name: str, settings: Any) -> list[Path]:
    paths = [directory / f"{name}.pdf", directory / f"{name}.png"]
    # Do not crop: preserve the configured physical text width in the PDF.
    fig.savefig(paths[0], metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(paths[1], dpi=_get(settings, "dpi"))
    plt.close(fig)
    return paths


def _model_legend(models: list[str]) -> list:
    return [
        Line2D(
            [],
            [],
            color=_style(model)[1],
            marker=_style(model)[2],
            linestyle=_style(model)[3],
            label=_style(model)[0],
            markersize=8,
        )
        for model in models
    ]


def _range(ax, position, low, high, color, *, horizontal=False, linewidth=5, alpha=0.25):
    if not (_number(low) and _number(high)):
        return
    if horizontal:
        ax.hlines(position, low, high, color=color, linewidth=linewidth, alpha=alpha)
    else:
        ax.vlines(position, low, high, color=color, linewidth=linewidth, alpha=alpha)


def _draw_summary(ax_coverage, ax_width, position, summary, model, *, horizontal=False):
    if not _visible(summary):
        return
    _label, color, marker, _line = _style(model)
    _range(
        ax_coverage,
        position,
        summary.get("observable_coverage_lower"),
        summary.get("observable_coverage_upper"),
        color,
        horizontal=horizontal,
    )
    ipcw = summary.get("ipcw_coverage")
    if _number(ipcw):
        ci = _ci(summary, "ipcw_coverage", "ipcw_coverage_ci")
        if ci is not None:
            _range(ax_coverage, position, *ci, color, horizontal=horizontal, linewidth=1.5, alpha=1)
        xy = (ipcw, position) if horizontal else (position, ipcw)
        ax_coverage.plot(
            *xy,
            marker=marker,
            linestyle="none",
            color=color,
            markersize=8,
            markeredgecolor="white",
            markeredgewidth=0.8,
        )
        if ci is not None and _partial_ci(summary, "ipcw_coverage"):
            ax_coverage.annotate(
                "a", xy, xytext=(4, 4), textcoords="offset points", fontsize=8, color=color
            )
    else:
        # Open endpoints show the bounds when IPCW cannot be estimated.
        for bound in ("observable_coverage_lower", "observable_coverage_upper"):
            if _number(summary.get(bound)):
                xy = (summary[bound], position) if horizontal else (position, summary[bound])
                ax_coverage.plot(
                    *xy,
                    marker=marker,
                    linestyle="none",
                    color=color,
                    markerfacecolor="none",
                    markersize=8,
                )
    width = summary.get("width_median_months")
    if _number(width):
        _range(
            ax_width,
            position,
            summary.get("width_q25_months"),
            summary.get("width_q75_months"),
            color,
            horizontal=horizontal,
        )
        ci = _ci(summary, "width_median_months", "width_median_ci_months")
        if ci is not None:
            _range(ax_width, position, *ci, color, horizontal=horizontal, linewidth=1.5, alpha=1)
        xy = (width, position) if horizontal else (position, width)
        ax_width.plot(
            *xy,
            marker=marker,
            linestyle="none",
            color=color,
            markersize=8,
            markeredgecolor="white",
            markeredgewidth=0.8,
        )
    elif summary.get("width_median_is_infinite"):
        # No fabricated finite marker: the annotation carries the infinite value.
        transform = ax_width.get_yaxis_transform() if horizontal else ax_width.get_xaxis_transform()
        xy = (0.98, position) if horizontal else (position, 0.98)
        ax_width.text(
            *xy,
            r"$\infty$",
            transform=transform,
            color=color,
            ha="right" if horizontal else "center",
            va="top",
        )


def _marginal_figure(
    rows: list[dict], context: tuple, directory: Path, settings: Any
) -> list[Path]:
    fig, axes = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(_get(settings, "width_inches"), _get(settings, "height_inches") * 1.55),
    )
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.1, top=0.82, hspace=0.17)
    models = sorted({row["model"] for row in rows})
    nominals = sorted({1 - row["marginal"]["alpha"] for row in rows})
    for j, model in enumerate(models):
        subset = sorted(
            (row for row in rows if row["model"] == model),
            key=lambda row: row["marginal"]["nominal_coverage"],
        )
        x = np.array([row["marginal"]["nominal_coverage"] for row in subset])
        # Slight deterministic offsets separate coincident intervals; tick labels
        # and the underlying nominal values always retain their exact levels.
        offset = (j - (len(models) - 1) / 2) * 0.004
        for nominal, row in zip(x, subset, strict=True):
            _draw_summary(*axes, nominal + offset, row["marginal"], model)
        width = [row["marginal"].get("width_median_months") for row in subset]
        if all(_number(value) for value in width):
            axes[1].plot(
                x + offset,
                width,
                color=_style(model)[1],
                linestyle=_style(model)[3],
                linewidth=_get(settings, "line_width"),
                zorder=1,
            )
    axes[0].plot([0, 1], [0, 1], color=REFERENCE, linewidth=1, label="Nominal")
    axes[0].set(ylabel="Coverage", ylim=(-0.02, 1.04))
    axes[1].set(
        ylabel="Interval width (months)", xlabel="Nominal coverage", ylim=(0, context[3] * 1.04)
    )
    axes[1].set_xticks(nominals, [f"{value:.0%}" for value in nominals])
    axes[1].set_xlim(min(nominals) - 0.025, max(nominals) + 0.025)
    fig.legend(
        handles=_model_legend(models),
        loc="upper center",
        ncol=1,
        frameon=False,
        bbox_to_anchor=(0.55, 0.995),
    )
    fig.text(
        0.16,
        0.84,
        "Coverage: thick range = observable bounds; symbol / whisker = IPCW / CI.\n"
        "Width: symbol = median; thick range = IQR; whisker = median CI.",
        fontsize=max(8, _get(settings, "font_size") - 1),
    )
    for ax in axes:
        ax.grid(axis="y", color="#dddddd", linewidth=0.5)
        ax.set_axisbelow(True)
    return _save(fig, directory, f"rq5_coverage_vs_nominal__{_context_name(context)}", settings)


def _subgroup_figure(
    rows: list[dict], context: tuple, column: str, directory: Path, settings: Any
) -> list[Path]:
    models = sorted({row["model"] for row in rows})
    alphas = sorted({row["marginal"]["alpha"] for row in rows}, reverse=True)
    levels = sorted({level for row in rows for level in row.get("subgroups", {}).get(column, {})})
    if not levels:
        return []
    height = max(_get(settings, "height_inches") * 1.55, 1.8 + len(levels) * 0.72)
    fig, axes = plt.subplots(
        2, len(alphas), squeeze=False, sharey=True, figsize=(_get(settings, "width_inches"), height)
    )
    fig.subplots_adjust(left=0.24, right=0.98, bottom=0.08, top=0.80, wspace=0.14, hspace=0.22)
    for col, alpha in enumerate(alphas):
        for j, model in enumerate(models):
            matching = [
                row for row in rows if row["model"] == model and row["marginal"]["alpha"] == alpha
            ]
            if not matching:
                continue
            groups = matching[0].get("subgroups", {}).get(column, {})
            for k, level in enumerate(levels):
                summary = groups.get(level, {"n_test": 0})
                position = k + (j - (len(models) - 1) / 2) * 0.20
                _draw_summary(axes[0, col], axes[1, col], position, summary, model, horizontal=True)
        axes[0, col].axvline(1 - alpha, color=REFERENCE, linewidth=1)
        axes[0, col].set(xlim=(-0.02, 1.04), xlabel=f"Coverage ({1 - alpha:.0%} nominal)")
        axes[1, col].set(xlim=(0, context[3] * 1.04), xlabel="Width (months)")
        axes[0, col].set_xticks([0, 0.5, 1])
        axes[1, col].set_xticks([0, context[3] / 2, context[3]])
        for ax in axes[:, col]:
            ax.set_yticks(range(len(levels)))
            ax.set_ylim(len(levels) - 0.5, -0.5)
            ax.grid(axis="x", color="#dddddd", linewidth=0.5)
            ax.set_axisbelow(True)
    first = next(row for row in rows if column in row.get("subgroups", {}))
    labels = []
    for level in levels:
        group = first["subgroups"][column].get(level, {"n_test": 0})
        suffix = "; suppressed" if group.get("status") == "suppressed_small_group" else ""
        labels.append(f"{level}\n(n={group.get('n_test', 0)}{suffix})")
    axes[0, 0].set_yticklabels(labels)
    axes[1, 0].set_yticklabels(labels)
    fig.legend(
        handles=_model_legend(models),
        loc="upper center",
        ncol=1,
        frameon=False,
        bbox_to_anchor=(0.60, 0.995),
    )
    fig.text(
        0.04,
        0.84,
        "Coverage: thick range = observable bounds; symbol / whisker = IPCW / CI.\n"
        "Width: symbol = median; thick range = IQR; whisker = median CI.\n"
        "Blank rows: small group suppressed or no observations.",
        fontsize=max(8, _get(settings, "font_size") - 1),
    )
    return _save(
        fig,
        directory,
        f"rq5_subgroup_coverage_width__{_context_name(context)}__{_slug(column)}",
        settings,
    )


def _cell(value: Any, digits: int = 2, *, infinite=False) -> str:
    if infinite:
        return r"\multicolumn{1}{c}{$\infty$}"
    return f"{value:.{digits}f}" if _number(value) else r"\multicolumn{1}{c}{--}"


def _interval_cell(value: Any, digits=2) -> str:
    if value is None or len(value) != 2 or not all(_number(v) for v in value):
        return "--"
    return f"[{value[0]:.{digits}f}, {value[1]:.{digits}f}]"


def _ci_cell(summary: dict, statistic: str, key: str, digits=2) -> str:
    ci = _ci(summary, statistic, key)
    value = _interval_cell(ci, digits)
    return value + (
        r"\textsuperscript{a}" if ci is not None and _partial_ci(summary, statistic) else ""
    )


def _table(rows: list[dict], context: tuple, column: str | None, directory: Path) -> Path:
    table_rows = []
    for row in sorted(rows, key=lambda x: (-x["marginal"]["alpha"], x["model"])):
        groups = (
            row.get("subgroups", {}).get(column, {}) if column else {"Overall": row["marginal"]}
        )
        for level, summary in sorted(groups.items()):
            table_rows.append(
                (
                    level,
                    row["model"],
                    row["marginal"]["nominal_coverage"],
                    row["marginal"]["n_cal"],
                    summary,
                )
            )
    prefix = _context_name(context)
    name = f"rq5_{'subgroups' if column else 'marginal'}__{prefix}"
    if column:
        name += f"__{_slug(column)}"
    path = directory / f"{name}.tex"
    preamble = [
        "% Generated from stored conformal_coverage.json; no metrics recomputed.",
        "% Requires booktabs, longtable, siunitx. Observable bounds are not latent-coverage CIs.",
        "% CI uncertainty is conditional on fixed model/calibration/train censoring fits.",
        "% Retained CIs after missing bootstrap draws are success-conditioned and approximate.",
        "% Width quartiles describe interval dispersion; median CI describes sampling uncertainty.",
        r"\begingroup",
        r"\setlength{\tabcolsep}{3pt}",
        r"\small",
    ]
    # Separate panels keep all estimates and uncertainty at thesis text width.
    coverage = [
        r"\begin{longtable}{@{}ll S[table-format=2.0] r r "
        r"S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] l@{}}",
        r"\toprule",
        r"Group & Model & {Nom.\%} & $n_{cal}$ & $n$ & {Lower} & {Upper} & {IPCW} & CI \\",
        r"\midrule\endhead",
    ]
    widths = [
        r"\begin{longtable}{@{}ll S[table-format=2.0] r S[table-format=3.1] l l "
        r"S[table-format=1.2] S[table-format=1.2]@{}}",
        r"\toprule",
        r"Group & Model & {Nom.\%} & $n$ & {Median} & IQR & Median CI & {Full} & {$\infty$} \\",
        r"\midrule\endhead",
    ]
    bound_intervals = [
        r"\begin{longtable}{@{}ll S[table-format=2.0] r l l@{}}",
        r"\toprule",
        r"Group & Model & {Nom.\%} & $n$ & Lower-bound CI & Upper-bound CI \\",
        r"\midrule\endhead",
    ]
    for level, model, nominal, n_cal, summary in table_rows:
        visible = _visible(summary)
        safe = summary if visible else {}
        group_label = _tex(level) + (
            r"\textsuperscript{s}" if summary.get("status") == "suppressed_small_group" else ""
        )
        model_label = {"cox": "Cox", "weibull": "Weibull", "rsf": "RSF", "lognormal": "Log-normal"}[
            model
        ]
        shared = [group_label, model_label, f"{100 * nominal:.0f}"]
        bound_intervals.append(
            " & ".join(
                [
                    *shared,
                    str(summary.get("n_test", 0)),
                    _ci_cell(safe, "observable_coverage_lower", "observable_coverage_lower_ci"),
                    _ci_cell(safe, "observable_coverage_upper", "observable_coverage_upper_ci"),
                ]
            )
            + r" \\"
        )
        coverage.append(
            " & ".join(
                [
                    *shared,
                    str(n_cal),
                    str(summary.get("n_test", 0)),
                    _cell(safe.get("observable_coverage_lower")),
                    _cell(safe.get("observable_coverage_upper")),
                    _cell(safe.get("ipcw_coverage")),
                    _ci_cell(safe, "ipcw_coverage", "ipcw_coverage_ci"),
                ]
            )
            + r" \\"
        )
        widths.append(
            " & ".join(
                [
                    *shared,
                    str(summary.get("n_test", 0)),
                    _cell(
                        safe.get("width_median_months"),
                        1,
                        infinite=safe.get("width_median_is_infinite", False),
                    ),
                    _interval_cell([safe.get("width_q25_months"), safe.get("width_q75_months")], 1),
                    _ci_cell(safe, "width_median_months", "width_median_ci_months", 1),
                    _cell(safe.get("full_support_fraction")),
                    _cell(safe.get("infinite_interval_fraction")),
                ]
            )
            + r" \\"
        )
    ending = [r"\bottomrule", r"\end{longtable}"]
    path.write_text(
        "\n".join(
            [
                *preamble,
                *coverage,
                *ending,
                r"\medskip",
                *widths,
                *ending,
                r"\medskip",
                *bound_intervals,
                *ending,
                r"\noindent\textsuperscript{s}Small group: estimates suppressed. "
                r"Widths in months; Full and $\infty$ are proportions.",
                r"\par\noindent\textsuperscript{a}CI conditions on successful bootstrap draws "
                r"and is approximate.",
                r"\endgroup",
                "",
            ]
        )
    )
    return path


def render_coverage(report: dict, output_dir: Path, cfg: Any) -> list[Path]:
    """Render only supplied aggregate stage-05 JSON to figures/ and tables/.

    Each analysis context gets a marginal coverage/width figure. Subgroup figures
    use columns for nominal levels and rows for coverage and width, so no alpha
    is silently selected for presentation. Suppressed small groups remain blank.
    """
    settings, conformal = _get(cfg, "figures"), _get(cfg, "conformal")
    configured_alphas = _get(conformal, "alphas")
    configured_horizons = _get(conformal, "horizons_months")
    contexts = defaultdict(list)
    for row in report.values():
        summary = row["marginal"]
        if (
            summary["alpha"] not in configured_alphas
            or summary["horizon_months"] not in configured_horizons
        ):
            raise ValueError("Coverage metrics do not match the configured alpha/horizon protocol.")
        _style(row["model"])
        context = (
            row["endpoint"],
            row["feature_set"],
            row["registered_method"],
            summary["horizon_months"],
        )
        contexts[context].append(row)
    directory = Path(output_dir)
    figures, tables = directory / "figures", directory / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    paths = []
    with plt.rc_context(
        {
            "font.size": _get(settings, "font_size"),
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        for context, rows in sorted(contexts.items()):
            paths.extend(_marginal_figure(rows, context, figures, settings))
            paths.append(_table(rows, context, None, tables))
            for column in sorted({column for row in rows for column in row.get("subgroups", {})}):
                paths.extend(_subgroup_figure(rows, context, column, figures, settings))
                paths.append(_table(rows, context, column, tables))
            caption = figures / f"rq5_caption__{_context_name(context)}.txt"
            endpoint, features, method, horizon = context
            counts = "; ".join(
                f"{_style(row['model'])[0]}: n_cal={row['marginal']['n_cal']}, "
                f"n_test={row['marginal']['n_test']}"
                for row in rows
                if row["marginal"]["alpha"] == max(configured_alphas)
            )
            confidence_levels = sorted({row["marginal"]["confidence_level"] for row in rows})
            confidence_label = ", ".join(f"{level:.0%}" for level in confidence_levels)
            caption.write_text(
                f"{endpoint.upper()}; {features}; {method}; "
                f"restricted target min(T,{horizon:g} months).\n"
                f"{counts}. Calibration n=0 denotes native model predictive intervals.\n"
                f"Confidence level: {confidence_label}.\n"
                "Thick coverage ranges are observable lower/upper bounds, "
                "not confidence intervals for latent survival coverage. "
                "Filled symbols and thin whiskers show IPCW coverage and its bootstrap CI "
                "under independent-censoring assumptions. "
                "Open endpoints identify observable bounds when IPCW is unavailable. "
                "Grey references show nominal coverage.\n"
                "Width symbols show the median, thick ranges the IQR, and thin whiskers "
                "the bootstrap median CI. Bootstrap CIs condition on fitted models, calibration "
                "and the train censoring fit; incomplete resampling yields "
                "success-conditioned approximate CIs. Marginal nominal levels "
                "do not guarantee subgroup coverage. Small-group estimates are suppressed. "
                "The annotation a marks CIs conditioned on successful resamples. "
                "See accompanying tables and conformal_coverage.json for statuses.\n"
            )
            paths.append(caption)
    return paths
