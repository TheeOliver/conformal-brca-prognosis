---
description: How predictive performance and coverage are measured, and the time-grid trap behind both.
paths:
  - "src/brca/evaluation/**"
  - "scripts/05_evaluate.py"
---

# Evaluation metrics

- **The time grid is shared and interior.** Build it once from
  `config.evaluation.time_grid_percentiles` as percentiles of observed follow-up, and keep it
  strictly inside the follow-up range of **both** the training and the evaluation set. IPCW
  weights explode near the tails, so a grid point past the last training event turns an
  IBS into a meaningless number without raising anything.
- Report together, never in isolation: Uno's C-index, Brier score, Integrated Brier Score, and
  time-dependent AUC where a specific horizon matters.
- **Coverage and width are reported as a pair**, at every alpha in `config.conformal.alphas`.
  Coverage alone cannot distinguish a good interval from an uninformative one.
- **Conditional coverage is required, not optional.** Break coverage down by
  `config.conformal.coverage_subgroups`. Marginal coverage can sit exactly at 90% while a
  PAM50 subtype is covered at 60% — that is an equity failure and a headline result.
- Attach uncertainty to every metric: bootstrap over the test set
  (`config.evaluation.bootstrap_replicates`) and report the interval, not just the point.
- A metric is computed **once** into `outputs/metrics/*.json`; figures, tables and prose read
  that file. Never recompute for a plot — two versions of one number is how a defence fails.
- Discrimination and calibration are different claims. A model may win on C-index and lose on
  Brier; say so rather than picking whichever is favourable. That contrast is RQ6.

Definitions and worked interpretation: `docs/survival-metrics.md`.
