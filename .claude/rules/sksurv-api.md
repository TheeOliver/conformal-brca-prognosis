---
description: scikit-survival API choices and the traps that silently produce wrong numbers.
paths:
  - "src/brca/models/**"
  - "src/brca/evaluation/**"
---

# scikit-survival API

Each of these is a real failure mode, not a style preference:

- **Use `sksurv.preprocessing.OneHotEncoder`, not sklearn's.** The sksurv one takes and returns
  a DataFrame with readable column names; sklearn's returns a bare array and destroys the
  column mapping needed to interpret Cox coefficients.
- **Use `concordance_index_ipcw`, not `concordance_index_censored`.** Harrell's C is biased
  upward when test censoring is heavy, which it is here. Uno's C needs the *training* labels
  for the censoring distribution:
  `concordance_index_ipcw(y_train, y_test, scores)`. Passing `y_test` in both positions is a
  silent wrong answer, not an error.
- For cross-validation inside train, wrap estimators with `as_concordance_index_ipcw_scorer`,
  `as_cumulative_dynamic_auc_scorer` or `as_integrated_brier_score_scorer` — plain sklearn
  scorers cannot read a structured survival label.
- `predict()` returns a **risk score** (higher = worse), not a time.
  `predict_survival_function()` returns `StepFunction` objects; evaluate them on the shared
  time grid rather than comparing them at arbitrary points.
- Keep the sksurv encoder inside the `make_pipeline`, so it is fit on train with everything
  else and cannot leak.

Which metric answers which research question: `docs/survival-metrics.md`.
