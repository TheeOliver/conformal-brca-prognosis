# Survival metrics

> Pulled in by: `run-experiment`, `compare-feature-sets` skills;
> `.claude/rules/evaluation-metrics.md`, `.claude/rules/sksurv-api.md`.

## Which metric answers which question

| Metric | Measures | Answers |
| --- | --- | --- |
| Uno's C-index (IPCW) | discrimination — ranking | RQ1, RQ3 |
| Brier score at t | calibration + discrimination at one horizon | RQ1, RQ6 |
| Integrated Brier Score | the above across the grid | RQ3, RQ6 |
| time-dependent AUC | discrimination at a clinical horizon | RQ1 |
| empirical coverage | validity of the interval | RQ5 |
| interval width | efficiency of the interval | RQ5, RQ6 |

Discrimination and calibration are **different claims**. A model can rank patients well while
its predicted probabilities are badly miscalibrated — good C-index, poor Brier. RQ6 exists
because of this dissociation, so never collapse them into "performance".

## Uno's C, not Harrell's

```python
from sksurv.metrics import concordance_index_ipcw

cindex = concordance_index_ipcw(y_train, y_test, risk_scores, tau=tau)[0]
```

Harrell's `concordance_index_censored` is biased upward when test censoring is heavy, which it
is in METABRIC. Uno's version reweights by the censoring distribution, which it estimates from
the **first** argument. Passing `y_test` in both positions estimates censoring from the wrong
sample and returns a plausible, wrong number with no warning.

## The time grid — the quiet failure

Build it once, from percentiles of observed follow-up, and keep it strictly **interior**:

```python
times = np.percentile(y_train["time"], config.evaluation.time_grid_percentiles)
```

IPCW weights are `1 / G(t)`, with `G` the censoring survival function. As `t` approaches the end
of follow-up, `G(t)` approaches zero and the weights explode. A grid point past the last
training event produces a meaningless IBS and raises nothing.

The grid must lie inside the follow-up range of **both** the training and the evaluation set,
and the **same grid** must serve every model and both feature sets — otherwise the comparison
table is incoherent.

## Uncertainty on every metric

Bootstrap over test patients (`config.evaluation.bootstrap_replicates`). When comparing the two
feature-set arms, bootstrap **paired** on the same patients: the arms are not independent
samples, and an unpaired comparison overstates the uncertainty of the difference.

A difference within bootstrap noise is not a difference. Say so.

## Reference

`references/notebooks_reference/Evaluating Survival Models guide (...).md` — quote the filename
in shell commands, and per `.claude/rules/notebooks.md` treat its code as lossy: read for the
API surface, never copy.
