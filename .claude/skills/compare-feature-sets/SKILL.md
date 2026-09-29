---
name: compare-feature-sets
description: Compare the clinical model against the clinical+PAM50 model on identical splits across prediction, calibration and uncertainty. Use when asked whether PAM50 or molecular information adds value, whether the molecular model is better, or anything touching RQ2.
---

# Clinical vs clinical + PAM50

This answers RQ2 — *does the molecular subtype improve the prognosis?* The comparison is only
meaningful if everything except the feature columns is held fixed.

Read `docs/experimental-design.md` and `docs/survival-metrics.md` first.

## 1. Verify the comparison is fair

Before computing anything, assert:

- Both arms used the **same** `data/processed/splits.json` — same patient IDs in each of the
  three sets. Different splits make any difference uninterpretable.
- Same seed, same time grid, same alphas, same preprocessing apart from the PAM50 columns.
- Both arms exist for **every** model (Cox, Bayesian AFT, RSF). A missing cell is a gap in the
  answer, not something to quietly omit from the table.
- Patients with missing PAM50 are handled identically in both arms — dropping them only from
  the molecular arm changes the population and invalidates the comparison outright.

## 2. Compare on all three axes

PAM50 can help on one axis and not another; that dissociation is the interesting part.

| Axis | Metric | Question |
| --- | --- | --- |
| Discrimination | Uno's C-index, time-dependent AUC | does it rank patients better? |
| Calibration | Brier, IBS | are the predicted probabilities better? |
| Uncertainty | interval width at fixed coverage | does it make prognosis *sharper*? |

The third axis is the one a standard comparison misses. If both arms reach 90% coverage but
the molecular arm's intervals are narrower, PAM50 bought **precision** even when the C-index
barely moved. That is a genuine RQ2 finding and it only shows up if width is reported.

## 3. Quantify the difference, not just the direction

- Paired bootstrap over the **same** test patients — the two arms are not independent samples,
  so an unpaired comparison overstates the uncertainty of the difference.
- Report the difference with its interval, not two separate point estimates side by side.
- Do not run a significance test per metric per model and then report the smallest p-value.
  State up front which comparison is primary.

## 4. Interpret honestly

- A gain within bootstrap noise is **not** a gain. Say so.
- PAM50 is partly redundant with grade, ER and HER2, which are already in the clinical arm —
  a small increment is the expected result, not a failure of the analysis.
- Use associative language: PAM50 improves *predictive performance*, never *patient outcomes*
  (`.claude/rules/scientific-claims.md`).

## 5. Emit

`outputs/metrics/rq2_feature_set_comparison.json`, then a booktabs table and a paired
coverage/width figure via `export-thesis-figure`. Write the conclusion — including a null
one — into `docs/experimental-design.md` under the RQ2 heading.
