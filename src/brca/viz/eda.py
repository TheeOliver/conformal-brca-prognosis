"""Export EDA figures from a saved training-EDA aggregate JSON only."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PALETTE = {
    "LumA": ("#2a78d6", "-"),
    "LumB": ("#eb6834", "--"),
    "Her2": ("#1baf7a", "-."),
    "Basal": ("#eda100", ":"),
    "Normal": ("#e87ba4", (0, (6, 2))),
    "claudin-low": ("#6b4c9a", (0, (3, 1, 1, 1))),
    "unknown": ("#8c564b", (0, (1, 1))),
    "positive": ("#2a78d6", "-"),
    "negative": ("#eb6834", "--"),
    "overall": ("#2a78d6", "-"),
}


def save_figure(fig, directory: Path, name: str, config) -> list[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    paths = [directory / f"{name}.pdf", directory / f"{name}.png"]
    fig.savefig(paths[0], metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(paths[1], dpi=config.figures.dpi)
    plt.close(fig)
    return paths


def export_training_eda(report: dict, directory: Path, config) -> list[Path]:
    plt.rcParams.update(
        {
            "font.size": config.figures.font_size,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    paths = []
    for name, summary in report["numeric"].items():
        fig, ax = plt.subplots(figsize=(config.figures.width_inches, config.figures.height_inches))
        edges = np.array(summary["histogram_edges"])
        ax.bar(
            edges[:-1],
            summary["histogram_counts"],
            width=np.diff(edges),
            align="edge",
            color="#2a78d6",
            edgecolor="white",
        )
        units = {
            "age_at_diagnosis": "Age at diagnosis (years)",
            "tumor_size": "Tumour size (mm)",
            "lymph_nodes_positive": "Positive lymph nodes (count)",
        }
        ax.set(xlabel=units[name], ylabel="Training patients")
        fig.subplots_adjust(left=0.14, right=0.97, bottom=0.18, top=0.97)
        paths.extend(save_figure(fig, directory, f"training_{name}_distribution", config))
    for endpoint, summaries in report["survival"].items():
        groups = {"overall": {"overall": summaries["overall"]}, **summaries["groups"]}
        for variable, curves in groups.items():
            n_curves = sum("survival" in curve for curve in curves.values())
            fig, (ax, risk_ax) = plt.subplots(
                2,
                1,
                sharex=True,
                figsize=(
                    config.figures.width_inches,
                    config.figures.height_inches + n_curves * 0.15,
                ),
                gridspec_kw={"height_ratios": [4, max(1, n_curves * 0.3)]},
            )
            fig.subplots_adjust(left=0.20, right=0.97, bottom=0.06, top=0.97, hspace=0.42)
            risk_rows, labels = [], []
            for level, curve in curves.items():
                if "survival" not in curve:
                    continue
                colour, linestyle = PALETTE.get(level, PALETTE["unknown"])
                time = np.r_[0, curve["times_months"]]
                survival = np.r_[1, curve["survival"]]
                ax.step(
                    time,
                    survival,
                    where="post",
                    label=f"{level} (n={curve['n']})",
                    color=colour,
                    linestyle=linestyle,
                    linewidth=config.figures.line_width,
                )
                ax.fill_between(
                    time,
                    np.r_[1, curve["lower"]],
                    np.r_[1, curve["upper"]],
                    step="post",
                    color=colour,
                    alpha=0.15 if variable == "overall" else 0.06,
                )
                censor_times = np.array(curve["censor_times_months"])
                if len(censor_times):
                    indices = np.searchsorted(time, censor_times, side="right") - 1
                    ax.plot(
                        censor_times,
                        survival[indices],
                        "|",
                        color="#767676",
                        markersize=3,
                        alpha=0.45,
                    )
                risk_rows.append(
                    [curve["at_risk"][str(t)] for t in config.eda.followup_horizons_months]
                )
                labels.append(level)
            ax.set(
                xlabel="Time (months)",
                ylabel=f"{endpoint.replace('_', ' ').capitalize()} survival",
                ylim=(0, 1.02),
            )
            ax.legend(loc="best", fontsize=config.figures.font_size - 1)
            ax.tick_params(labelbottom=True)
            risk_ax.set_ylim(0, len(labels) + 1.4)
            risk_ax.axis("off")
            for j, label in enumerate(["At risk", *labels]):
                y = len(labels) + 0.7 - j
                risk_ax.text(
                    -0.02,
                    y,
                    label,
                    transform=risk_ax.get_yaxis_transform(),
                    ha="right",
                    va="center",
                    color="#767676",
                    fontsize=config.figures.font_size - 1,
                )
                values = config.eda.followup_horizons_months if j == 0 else risk_rows[j - 1]
                for t, value in zip(config.eda.followup_horizons_months, values, strict=True):
                    risk_ax.text(
                        t,
                        y,
                        str(value),
                        ha="center",
                        va="center",
                        color="#767676",
                        fontsize=config.figures.font_size - 1,
                    )
            paths.extend(save_figure(fig, directory, f"training_{endpoint}_{variable}_km", config))
    return paths
