---
description: The three-way train/calibration/test protocol and the leakage rules that protect it.
paths:
  - "src/brca/**"
  - "scripts/**"
---

# Splits & leakage

The split is drawn **once**, by `scripts/02_make_splits.py`, from `config.seed`, and frozen in
`data/processed/splits.json` as patient IDs. Everything downstream reads that file.

- **train** fits models. **calibration** produces conformal quantiles and nothing else.
  **test** is read once, in `scripts/05_evaluate.py`, and never used to choose anything.
- Never call `train_test_split` anywhere else. Never re-split "just to check" — it silently
  invalidates every cached fit and every calibrated interval.
- Every fitted transform (imputer, scaler, encoder, PAM50 level set, any variable selection) is
  `fit` on **train only**, then `transform`ed onto calibration and test. Fitting on the
  concatenation is the leak that quietly inflates every number downstream.
- Categorical levels come from train. An unseen PAM50 level at test time maps to an explicit
  `unknown` column; it must not add a column and change the design-matrix width.
- Conformal validity rests on calibration and test being **exchangeable**. Anything that treats
  them differently — different imputation, filtering, or time origin — breaks the coverage
  guarantee itself, not merely the accuracy.
- Never tune a hyperparameter, threshold or time grid on calibration or test. If you need model
  selection, cross-validate **inside train**.

Protocol and rationale: `docs/experimental-design.md`.
