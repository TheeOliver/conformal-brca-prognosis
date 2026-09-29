---
name: add-survival-model
description: Add a new survival model (or a variant of an existing one) behind the shared interface so it drops into the pipeline, the conformal layer and the evaluation table unchanged. Use when asked to add, implement or wire up a model such as a penalised Cox, an AFT variant, a gradient-boosted survival model, or a new prior specification.
---

# Add a survival model

Every model in this project is interchangeable from the pipeline's point of view. Honour the
interface and nothing downstream needs to change.

## 1. Read the contract first

`docs/model-interface.md` defines what every model must provide:

- `fit(X_train, y_train) -> self`
- `predict_risk(X) -> np.ndarray` (higher = worse; used for the C-index)
- `predict_survival_function(X, times) -> np.ndarray` shape `(n_samples, n_times)`
- `save(path)` / `load(path)`

`predict_survival_function` is the one the conformal layer and the Brier score both consume,
so it must be correct on the **shared time grid**, not on the model's own internal grid.

## 2. Implement

1. New module under `src/brca/models/<name>.py`. One model per file.
2. Read the rules that apply: `.claude/rules/sksurv-api.md` for anything scikit-survival,
   `.claude/rules/bayesian-pymc.md` and `.claude/rules/mcmc-diagnostics.md` for anything PyMC.
3. All hyperparameters come from `config/default.yaml` under a key named after the model.
   No literals — see `.claude/rules/reproducibility.md`.
4. Accept an injected `rng`; never touch global numpy random state.
5. If the model is Bayesian, it must also expose the diagnostics gate and persist
   InferenceData — a point estimate from a Bayesian model is not the deliverable.

## 3. Test before wiring

Add `tests/test_models_<name>.py` covering:

- fits on `tests/fixtures/synthetic_metabric.csv` without error;
- `predict_survival_function` is **monotone non-increasing** in time and stays within `[0, 1]`
  — the most common silent bug in a new survival model;
- risk ordering is sane: a synthetic high-risk patient outranks a low-risk one;
- two fits with the same seed give identical predictions.

Mark anything that samples as `@pytest.mark.slow` so the Stop gate stays fast.

## 4. Wire in

1. Register in the model registry in `src/brca/models/__init__.py`.
2. Add its config block.
3. Run `make data split` then the new model's stage on synthetic data.
4. Confirm it appears in the evaluation table for **both** feature sets — a model added to
   only the clinical arm breaks the RQ2 comparison.

## 5. Document

Add a row to the model table in `docs/experimental-design.md`: what it assumes, what it buys,
and how it differs from the three baseline models. Then run the `leakage-auditor` agent.
