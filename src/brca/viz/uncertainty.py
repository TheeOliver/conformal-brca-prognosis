"""Draw hypothetical-profile Bayesian uncertainty only from its saved metrics JSON."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def export_hypothetical_uncertainty(
    metrics_path: Path, output_directory: Path, config
) -> list[Path]:
    """Export a two-panel PDF/PNG per hypothetical profile, without recomputing numbers."""
    report = json.loads(Path(metrics_path).read_text())
    if report.get("patient_records_included") is not False:
        raise ValueError("The uncertainty figure requires explicitly hypothetical profiles")
    output_directory.mkdir(parents=True, exist_ok=True)
    styles = {
        "clinical": ("Clinical", "#2a78d6", "-", "o"),
        "clinical_molecular": ("Clinical + PAM50", "#eb6834", "--", "s"),
    }
    times = report["time_grid_months"]
    artifacts = []
    with plt.rc_context({"font.size": config.figures.font_size, "pdf.fonttype": 42}):
        for profile in report["profiles"]:
            if not profile.get("hypothetical"):
                raise ValueError("A figure profile lacks its hypothetical designation")
            fig, axes = plt.subplots(
                1,
                2,
                figsize=(config.figures.width_inches, config.figures.height_inches),
                layout="constrained",
            )
            for feature, (label, color, line, marker) in styles.items():
                values = profile["models"][feature]
                axes[0].plot(
                    times,
                    values["survival_mean"],
                    label=label,
                    color=color,
                    linestyle=line,
                    linewidth=config.figures.line_width,
                )
                axes[0].fill_between(
                    times,
                    values["survival_credible_lower"],
                    values["survival_credible_upper"],
                    color=color,
                    alpha=0.15,
                )
                offset = -0.12 if feature == "clinical" else 0.12
                for j, interval in enumerate(values["predictive_intervals"]):
                    y = j + offset
                    axes[1].plot(
                        [interval["lower_months"], interval["upper_months"]],
                        [y, y],
                        color=color,
                        linestyle=line,
                        linewidth=config.figures.line_width,
                    )
                    axes[1].plot(
                        values["predictive_median_months"],
                        y,
                        color=color,
                        marker=marker,
                        markersize=8,
                        linestyle="none",
                    )
            axes[0].set(
                xlabel="Follow-up (months)", ylabel="Predicted survival probability", ylim=(0, 1)
            )
            axes[0].legend(loc="lower left", frameon=False)
            intervals = profile["models"]["clinical"]["predictive_intervals"]
            axes[1].set_yticks(
                np.arange(len(intervals)),
                [f"{100 * row['nominal_probability']:g}%" for row in intervals],
            )
            axes[1].set(
                xlabel="Predicted survival time (months)", ylabel="Central predictive probability"
            )
            axes[1].set_xlim(left=0)
            axes[1].margins(y=0.3)
            for ax in axes:
                ax.spines[["top", "right"]].set_visible(False)
                ax.tick_params(labelsize=config.figures.font_size)
            stem = output_directory / f"rq4_hypothetical_uncertainty_{profile['profile_id']}"
            for extension in ("pdf", "png"):
                path = stem.with_suffix(f".{extension}")
                metadata = {"Title": ""}
                if extension == "pdf":
                    metadata.update({"CreationDate": None, "ModDate": None})
                fig.savefig(path, dpi=config.figures.dpi, metadata=metadata)
                artifacts.append(path)
            plt.close(fig)
    return artifacts
