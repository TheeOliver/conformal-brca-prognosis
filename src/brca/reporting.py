"""Traceable research-question answers assembled from saved aggregate results."""

from __future__ import annotations

import json
from pathlib import Path

from brca.pipeline import write_json


def _number(value, digits=3):
    return "unavailable" if value is None else f"{value:.{digits}f}"


def _estimate(value, digits=3):
    interval = f"[{_number(value['ci_lower'], digits)}, {_number(value['ci_upper'], digits)}]"
    suffix = "" if value["status"] == "ok" else f" ({value['status']})"
    return f"{_number(value['estimate'], digits)} {interval}{suffix}"


def _difference_interpretation(summary, lower_better=False):
    lo, hi = summary["ci_lower"], summary["ci_upper"]
    if lo is None or hi is None:
        return "Uncertainty interval unavailable; no directional conclusion."
    if lo <= 0 <= hi:
        return "The interval includes zero; the direction of added predictive value is uncertain."
    favors = hi < 0 if lower_better else lo > 0
    label = "clinical + PAM50" if favors else "clinical alone"
    return (
        f"The paired interval favors {label}; "
        "comparisons are exploratory without multiplicity adjustment."
    )


def write_research_report(outputs: Path, config) -> list[Path]:
    metrics = outputs / "metrics"
    performance = json.loads((metrics / "performance.json").read_text())
    coverage = json.loads((metrics / "conformal_coverage.json").read_text())
    comparisons = json.loads((metrics / "rq2_feature_set_comparison.json").read_text())
    fits = json.loads((outputs / "models/fit_index.json").read_text())
    sensitivity = json.loads((metrics / "bayesian_sensitivity.json").read_text())
    primary = config.outcome.primary_endpoint
    p = performance[primary]
    first_model = next(iter(p["models"].values()))
    grid_start = first_model["brier"][0]["time_months"]
    text = [
        "# METABRIC thesis results",
        "",
        "This report is generated from saved aggregate JSON; it does not recompute metrics.",
        f"Primary endpoint: {primary}. OS is a sensitivity analysis. "
        f"Training n={p['n_training']}; test n={p['n_test']}. "
        f"The predictive metric grid spans {grid_start:.3f} to "
        f"{p['evaluation_horizon_months']:g} months; IBS integrates over "
        "that interval, not from zero.",
        "DSS describes net cause-specific survival under censoring assumptions, not crude "
        "breast-cancer mortality. All analyses are prognostic.",
        "",
        "## RQ1 — Clinical predictive performance",
        "",
        "| Model | Uno C (95% CI) | IBS (95% CI) | Mean dynamic AUC (95% CI) |",
        "| --- | --- | --- | --- |",
    ]
    for key, row in p["models"].items():
        if key.endswith("__clinical"):
            text.append(
                f"| {key.removesuffix('__clinical')} | {_estimate(row['uno_c'])} | "
                f"{_estimate(row['ibs'])} | {_estimate(row['mean_auc'])} |"
            )
    text += [
        "",
        "Evidence: `outputs/metrics/performance.json`.",
        "",
        "## RQ2 — Added predictive value of PAM50",
        "",
        "Differences are clinical + PAM50 minus clinical on identical patients.",
        "",
        "| Model | Paired Δ Uno C (95% CI) | Paired Δ IBS (95% CI) |",
        "| --- | --- | --- |",
    ]
    for model, row in comparisons["performance"][primary].items():
        text.append(f"| {model} | {_estimate(row['uno_c'])} | {_estimate(row['ibs'], 5)} |")
    for model, row in comparisons["performance"][primary].items():
        text += ["", f"{model}, primary Uno C: {_difference_interpretation(row['uno_c'])}"]
    text += [
        "",
        "Width and coverage changes also have paired intervals in "
        "`outputs/metrics/rq2_feature_set_comparison.json`; predictive gains alone do "
        "not establish better uncertainty.",
        "",
        "## RQ3 — Classical, Bayesian and forest models",
        "",
        "| Model / features | Uno C (95% CI) | IBS (95% CI) |",
        "| --- | --- | --- |",
    ]
    ranking = sorted(
        p["models"], key=lambda key: p["models"][key]["uno_c"]["estimate"], reverse=True
    )
    for key in ranking:
        row = p["models"][key]
        text.append(f"| {key} | {_estimate(row['uno_c'])} | {_estimate(row['ibs'])} |")
    text += [
        "",
        f"The largest primary Uno C point estimate belongs to {ranking[0]}. "
        "This descriptive ordering is not a test of superiority between model "
        "families. RSF was deliberately untuned.",
    ]
    ph_diagnostics = {}
    bayesian_diagnostics = {}
    for key, record in (fits["models"] | fits["sensitivities"]).items():
        report = json.loads(Path(record["metrics_path"]).read_text())
        if record["model"] == "cox":
            ph_diagnostics[key] = report["proportional_hazards"]
        elif record["model"] == "weibull":
            bayesian_diagnostics[key] = report
    text += ["", "Cox proportional-hazards score checks (training only, unadjusted):"]
    for key, diagnostic in ph_diagnostics.items():
        global_test = diagnostic.get("global", {})
        p_value = global_test.get("p_value")
        formatted_p = "unavailable" if p_value is None else f"{p_value:.3g}"
        text.append(
            f"- {key}: global p={formatted_p}; flagged={global_test.get('flagged', 'unavailable')}."
        )
    flagged = sum(bool(d.get("global", {}).get("flagged")) for d in ph_diagnostics.values())
    text += [
        "",
        f"Global PH departures were detected in {flagged} of {len(ph_diagnostics)} "
        "Cox fits. They remain prognostic benchmarks, with this assumption limitation explicit. "
        "Continuous coefficients are per training standard deviation, "
        "not per raw measurement unit.",
    ]
    text += [
        "",
        "A flagged PH check limits the proportional-hazards interpretation. "
        "Weibull AFT also imposes proportional hazards; adding AFT terminology does "
        "not remove that restriction.",
        "",
        "## RQ4 — Bayesian uncertainty",
        "",
        "Survival credible bands describe posterior uncertainty in S(t|x). Individual event-time "
        "predictive intervals additionally include variation of event times under the "
        "fitted model. "
        "Neither is a conformal coverage guarantee. Hypothetical profiles are "
        "constructed from training "
        "summaries and are not patient records.",
        "Holding modal receptors fixed while changing subtype can create rare joint predictor "
        "combinations. A subtype's training count alone does not establish support for that "
        "hypothetical combination; these illustrations are not validated individual forecasts.",
        "Unrestricted predictive tails can extend beyond plausible human lifetimes, especially "
        "for net DSS. Those are consequences of parametric extrapolation with no cure fraction "
        "or competing-risk distribution; they must not be read as literal clinical forecasts.",
        "",
        "| Fit | Max R-hat | Min bulk ESS | Min tail ESS | Divergences | Max training "
        "PPC discrepancy |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for key, report in bayesian_diagnostics.items():
        d, workflow = report["diagnostics"], report["workflow"]
        text.append(
            f"| {key} | {d['max_rhat']:.4f} | {d['min_ess_bulk']:.0f} | "
            f"{d['min_ess_tail']:.0f} | {d['divergences']} | "
            f"{max(workflow['training_ppc']['absolute_survival_discrepancy']):.3f} |"
        )
    text += [
        "",
        "Sensitivity: maximum absolute differences in fitted survival across training "
        "profiles and registered times:",
    ]
    for key, row in sensitivity.items():
        text.append(f"- {key}: {row['maximum_absolute_survival_difference']:.3f}.")
    text += [
        "",
        "Inspect the saved training PPC overlays: a passing MCMC gate establishes "
        "sampling adequacy, "
        "not a correct survival family. The alternative lognormal family and wider "
        "prior were checked "
        "on training patients only. Evidence: `bayes_diagnostics_*.json`, "
        "`bayesian_sensitivity.json` "
        "and the hypothetical-profile report in `outputs/metrics/`.",
        "",
        "## RQ5 — Coverage and efficiency under censoring",
        "",
        "Latent event times for censored patients are unknown. The observable lower/upper bounds "
        "enclose realised test coverage; they are not confidence limits for an identified coverage "
        "probability. IPCW coverage is a separate estimate requiring marginal "
        "independent censoring "
        "and adequate censoring support. The IPCW-calibrated method is an approximate sensitivity "
        "analysis with no asserted finite-sample guarantee.",
        "",
        "The table below shows primary DSS conservative intervals for min(T, horizon), at every "
        "registered nominal level. Coverage is always paired with width.",
        "",
        "| Model / features | Horizon (months) | Nominal | n calibration | Observable "
        "coverage bounds | IPCW coverage | Median width [IQR], months |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    primary_rows = []
    for key, row in coverage.items():
        if row["endpoint"] != primary or row["registered_method"] != "conservative":
            continue
        m = row["marginal"]
        primary_rows.append((key, row))
        text.append(
            f"| {row['model']} / {row['feature_set']} | {m['horizon_months']:g} | "
            f"{m['nominal_coverage']:.0%} | {m['n_cal']} | "
            f"[{m['observable_coverage_lower']:.3f}, "
            f"{m['observable_coverage_upper']:.3f}] | {_number(m['ipcw_coverage'])} | "
            f"{m['width_median_months']:.1f} [{m['width_q25_months']:.1f}, "
            f"{m['width_q75_months']:.1f}] |"
        )
    support = next(
        row["marginal"]
        for _, row in primary_rows
        if row["marginal"]["horizon_months"] == config.conformal.primary_horizon_months
    )
    if support["ipcw_status"] == "unavailable_censoring_support":
        censoring_survival = support["ipcw_diagnostics"]["censoring_survival_at_horizon"]
        text += [
            "",
            f"At {config.conformal.primary_horizon_months:g} months, estimated training "
            f"censoring survival is {_number(censoring_survival, 5)}, below the registered "
            f"support threshold {config.conformal.min_censoring_survival:g}. "
            "IPCW coverage is unavailable there. The IPCW-calibration sensitivity falls back "
            "to the entire restricted support, which is uninformative; it is not evidence of "
            "a useful calibrated interval. The 120-month analysis was prespecified.",
        ]
    text += [
        "",
        "Uncertainty intervals for each coverage bound, IPCW estimate and width statistic, "
        "plus PAM50/ER subgroup estimates, are in `conformal_coverage.json` and the "
        "exported tables. "
        "Groups smaller than the configured minimum are suppressed. Small remaining "
        "groups can still "
        "have wide uncertainty intervals. No subgroup coverage guarantee is inferred "
        "from marginal calibration.",
        "",
        "## RQ6 — Does the strongest discriminator give the most useful uncertainty?",
        "",
        f"Compare the largest Uno C point estimate ({ranking[0]}) with all interval widths above. "
        "A width close to the whole restricted support represents weak individual information even "
        "when coverage bounds are favorable. When bounds span nominal coverage or "
        "IPCW is unavailable, "
        "the data do not identify a uniquely best-calibrated model. Discrimination "
        "rankings therefore "
        "cannot by themselves establish which model has the most trustworthy uncertainty.",
        "",
        "## Overall-survival sensitivity",
        "",
        "| Model / features | Uno C (95% CI) | IBS (95% CI) |",
        "| --- | --- | --- |",
    ]
    for key, row in performance["overall"]["models"].items():
        text.append(f"| {key} | {_estimate(row['uno_c'])} | {_estimate(row['ibs'])} |")
    text += [
        "",
        "## Limits of inference",
        "",
        "One internal split, no external validation; historical treatment practice; "
        "missing-predictor "
        "imputation assumptions; possible dependent censoring and competing deaths; "
        "parametric family "
        "and PH misspecification; limited tail support; exploratory multiplicity. "
        "Test-patient bootstrap "
        "holds primary event-stratum counts fixed and conditions on fitted models, "
        "training censoring "
        "estimates and calibration. It omits model-fitting and split-selection "
        "variability. CIs with "
        "undefined replicates omitted are explicitly conditional and approximate; "
        "failed gates leave "
        "CIs unavailable. Full follow-up used for fitting must not be confused with "
        "restricted evaluation.",
        "",
        "Source/config archives, exact hashes and SLURM jobs are recorded in `outputs/manifests/`.",
    ]
    text += [
        "",
        "## RQ6 joint comparison at the primary restricted horizon",
        "",
        "| Model / features | Uno C (95% CI) | Nominal | Observable coverage bounds | "
        "Median width [95% CI], months | Full-support fraction |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in primary_rows:
        m = row["marginal"]
        if m["horizon_months"] != config.conformal.primary_horizon_months:
            continue
        model_key = f"{row['model']}__{row['feature_set']}"
        interval = m["width_median_ci_months"] or [None, None]
        text.append(
            f"| {model_key} | {_estimate(p['models'][model_key]['uno_c'])} | "
            f"{m['nominal_coverage']:.0%} | [{m['observable_coverage_lower']:.3f}, "
            f"{m['observable_coverage_upper']:.3f}] | {m['width_median_months']:.1f} "
            f"[{_number(interval[0], 1)}, {_number(interval[1], 1)}] | "
            f"{m['full_support_fraction']:.1%} |"
        )
    confidence = f"{config.evaluation.confidence_level:.0%}"
    text += [
        "",
        "## Subgroup signals requiring caution",
        "",
        "Below are primary-horizon subgroup cells whose observable upper coverage bound is "
        "below nominal on this test sample. Their bootstrap intervals concern the observable "
        "upper bound; intervals including nominal do not establish population undercoverage. "
        "These exploratory comparisons are not multiplicity-adjusted.",
        "",
        "| Model / features | Group | n | Nominal | Observable bounds | "
        "Upper-bound 95% CI | Median width, months |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in primary_rows:
        m = row["marginal"]
        if m["horizon_months"] != config.conformal.primary_horizon_months:
            continue
        for column, groups in row["subgroups"].items():
            for level, group in groups.items():
                upper = group.get("observable_coverage_upper")
                if upper is None or upper >= m["nominal_coverage"]:
                    continue
                ci = group.get("observable_coverage_upper_ci") or [None, None]
                text.append(
                    f"| {row['model']} / {row['feature_set']} | {column}: {level} | "
                    f"{group['n_test']} | {m['nominal_coverage']:.0%} | "
                    f"[{group['observable_coverage_lower']:.3f}, {upper:.3f}] | "
                    f"[{_number(ci[0])}, {_number(ci[1])}] | "
                    f"{group['width_median_months']:.1f} |"
                )
    text = [line.replace("95% CI", f"{confidence} CI") for line in text]
    report_path = outputs / "reports/thesis_results.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(text) + "\n")
    answers = {
        "primary_endpoint": primary,
        "descriptive_uno_c_ranking": ranking,
        "RQ1": {
            "evidence": "performance.json",
            "models": {k: v for k, v in p["models"].items() if k.endswith("__clinical")},
        },
        "RQ2": {
            "evidence": "rq2_feature_set_comparison.json",
            "paired_performance": comparisons["performance"][primary],
        },
        "RQ3": {"evidence": "performance.json", "ph_diagnostics": ph_diagnostics},
        "RQ4": {
            "evidence": ["bayes_diagnostics_*.json", "bayesian_sensitivity.json"],
            "all_diagnostic_gates_passed": all(
                x["diagnostics"]["passed"] for x in bayesian_diagnostics.values()
            ),
        },
        "RQ5": {"evidence": "conformal_coverage.json", "primary_conservative": dict(primary_rows)},
        "RQ6": {
            "evidence": ["performance.json", "conformal_coverage.json"],
            "latent_coverage_not_fully_observed": True,
        },
    }
    answer_path = metrics / "research_answers.json"
    write_json(answer_path, answers)
    return [report_path, answer_path]
