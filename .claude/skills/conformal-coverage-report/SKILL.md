---
name: conformal-coverage-report
description: Calibrate conformal prediction intervals at each alpha and report empirical coverage and interval width, marginally and by subgroup. Use when asked about coverage, whether intervals achieve 80/90/95%, validity versus efficiency, interval width, calibration of the conformal layer, or RQ5.
---

# Conformal coverage report

This answers RQ5 — *do the intervals actually achieve their declared coverage?* — and feeds
RQ6. Read `docs/conformal-methods.md` first; the rules in `.claude/rules/conformal.md` apply
throughout.

## 1. Check calibration integrity before computing anything

A coverage number from a contaminated calibration set is worse than no number.

1. Load `data/processed/splits.json`. Assert the calibration and test ID sets are **disjoint**
   and that neither overlaps train.
2. Confirm the underlying model was fit on **train only** — check the fit's manifest, not the
   code's intent.
3. Record `n_cal`. The guarantee is finite-sample; coverage reported without `n` is
   uninterpretable.

## 2. Compute the calibrated quantile, per alpha

For each `alpha` in `config.conformal.alphas` (0.20, 0.10, 0.05):

```
k = ceil((n_cal + 1) * (1 - alpha))
q = sorted(scores)[k - 1]          # if k > n_cal the interval is INFINITE
```

Not `np.quantile`. If `k > n_cal`, the calibration set is too small for that alpha — return an
infinite interval and say so; do not clip to the maximum score, which fabricates coverage.

## 3. Apply to test — once

Form intervals on the test set. This is the only time the test set is touched.

## 4. Measure both axes

For every alpha, report **as a pair**:

- **Validity** — empirical coverage, with a binomial confidence interval. Coverage of 0.88 on
  200 test patients is not distinguishable from 0.90; say that rather than calling it a miss.
- **Efficiency** — median and IQR of interval width, plus the proportion of infinite intervals.

## 5. Conditional coverage — required, not optional

Break coverage down by `config.conformal.coverage_subgroups` (PAM50 subtype, ER status).
Split conformal guarantees **marginal** coverage only: overall coverage can sit at exactly
90% while basal-like patients are covered at 60%. Finding that is a headline result about
trustworthiness and equity, not a footnote. Report subgroup `n` alongside — small subgroups
give noisy coverage and should not be over-read.

## 6. Contrast with Bayesian intervals

Put conformal intervals next to the Bayesian predictive intervals at the same nominal level.
The distinction is central to the thesis: the Bayesian interval says what the *model* believes
given its priors; the conformal interval is an empirical claim checked against held-out data.
Where they disagree, that gap *is* the result.

## 7. Emit

Write `outputs/metrics/conformal_coverage.json`, then build the table and figure via the
`export-thesis-figure` skill. Coverage figures always show the nominal level as a reference
line and pair coverage with width.

Finally, run the `conformal-validity-auditor` agent before any number leaves this repo.
