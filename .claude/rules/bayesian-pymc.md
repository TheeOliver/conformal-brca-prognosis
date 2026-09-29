---
description: PyMC idioms for the Weibull/AFT survival model - censoring, out-of-sample prediction, persistence.
paths:
  - "src/brca/models/bayesian*.py"
  - "src/brca/models/**"
---

# Bayesian model construction (PyMC)

Aliases: `import pymc as pm`, `import arviz as az`, `import pytensor.tensor as pt`.

- **Censoring uses `pm.Censored`** — one idiom, project-wide. Hand-rolling the log-survival
  term with `pm.Potential` is equivalent when correct and silently wrong when not. Reach for
  it only if `pm.Censored` genuinely cannot express the likelihood, and say why in a comment.
- Wrap predictors in `pm.Data(...)` so one fitted model can score calibration and test through
  `pm.set_data`. Rebuilding the model around new data re-estimates it and destroys the
  train-only guarantee in `.claude/rules/splits-and-leakage.md`.
- Sampler settings come from `config.mcmc` — `chains`, `tune`, `draws`, `target_accept`, and
  `random_seed` derived from `config.seed`. Never literals.
- **Priors are justified in `docs/bayesian-workflow.md`** and checked with a prior predictive
  draw before any posterior sampling. A prior that looks vague on an AFT log-scale is not vague
  on the time scale — it can put real mass on survival times of 10^6 months.
- Standardise continuous predictors (fit on train) before sampling, otherwise the sampler
  struggles and a shared prior means something different for each covariate.
- Persist the fit as InferenceData to `outputs/models/<name>.nc`. Never re-sample to regenerate
  a plot — load the stored fit.
- Sampling is compute: it runs under `sbatch`, never on the login node.

What makes a fit reportable: `.claude/rules/mcmc-diagnostics.md`.
