# Conformal prediction for survival data

> Pulled in by: `conformal-coverage-report` skill; `.claude/rules/conformal.md`;
> `conformal-validity-auditor` agent.

## Split (inductive) conformal prediction

1. Fit a model on **train**.
2. Compute a nonconformity score `s_i` for every **calibration** point.
3. Take the calibrated quantile `q` (below).
4. For a test point, the prediction set is everything whose score would be `<= q`.

Guarantee: `P(Y_test in C(X_test)) >= 1 - alpha`, over the draw of calibration and test
together. It is **marginal** — averaged over patients, not guaranteed for any subgroup.

## The finite-sample quantile — get this right

```python
k = math.ceil((n_cal + 1) * (1 - alpha))
q = np.sort(scores)[k - 1]  # k > n_cal  =>  interval is INFINITE
```

Not `np.quantile(scores, 1 - alpha)`. The `(n+1)` correction is what makes the guarantee
exact in finite samples; omitting it undercovers by about `1/n`. With `n_cal` around 500 that
is roughly 0.2% — undetectable by eye, and it means the stated guarantee is false.

When `k > n_cal` the calibration set is too small for that alpha. Return an infinite interval
and report the proportion. Clipping to `max(scores)` fabricates coverage.

## What censoring breaks

Plain CP assumes the outcome is observed. For a censored patient the true survival time is
unknown, so the nonconformity score is not computable. Approaches in the literature:

| Approach | Source | Idea |
| --- | --- | --- |
| Conformalised survival analysis | `candes_2023_conformalized_survival_analysis.md` | lower predictive bound on survival time; covariate-shift reweighting for censoring |
| Adaptive cut-offs | `gui_2024_conformalized_survival_adaptive_cutoffs.md` | data-driven threshold, better efficiency |
| IPCW-weighted scores | `qin_2024_conformal_predictive_intervals_survival.md` | reweight calibration by censoring probability |
| CONFIDE | `dang_2026_confide_conformal_causal_competing_risks.md` | competing risks; the one Python implementation |
| Beyond exchangeability | `barber_2023_conformal_prediction_beyond_exchangeability.md` | what survives when exchangeability fails |

**Implementation note.** `cfsurvival` (Candès) and Gui's adaptive cut-offs are **R**. Using
them means reimplementing, not importing. This is why `src/brca/conformal/` is hand-written and
why every score carries a coverage unit test.

## Validity vs efficiency

- **Validity** — does empirical coverage match nominal? Report with a binomial CI.
- **Efficiency** — how narrow? Median and IQR of width, plus the share of infinite intervals.

An interval spanning all follow-up trivially achieves 100% coverage and tells a clinician
nothing. Both axes, always, together.

## Conditional coverage

Marginal coverage can sit at exactly 90% while a PAM50 subtype is covered at 60%. Report
coverage by `config.conformal.coverage_subgroups` with subgroup `n`. Split CP does **not**
promise conditional coverage — measuring where it fails is a result, not a bug report.

## TODO

- [ ] Record which adaptation was implemented and what it assumes about censoring.
