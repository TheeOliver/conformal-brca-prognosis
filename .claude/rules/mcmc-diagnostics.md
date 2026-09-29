---
description: The acceptance gate every MCMC fit must pass before its numbers may be reported.
paths:
  - "src/brca/models/bayesian*.py"
  - "src/brca/evaluation/**"
  - "scripts/03_fit_models.py"
---

# MCMC diagnostics gate

A fit failing any check is **not reportable**. Thresholds live in `config.mcmc.diagnostics`;
assert them in code and **raise** — never warn and continue.

| Check | Threshold | Why it matters |
| --- | --- | --- |
| max R-hat | <= 1.01 | chains have not mixed |
| bulk ESS | >= 400 | posterior means are noise |
| tail ESS | >= 400 | **credible-interval endpoints are noise** |
| divergences | 0 | geometry unexplored, posterior biased |

Tail ESS carries unusual weight: this thesis reports *intervals*, and endpoints are tail
quantities. Required for every Bayesian fit, in order:

1. **Prior predictive check** — simulated times land on a clinically plausible scale.
2. Sample, then assert the table above.
3. **Posterior predictive check** — overlay posterior predictive survival on train Kaplan-Meier.
4. **Prior sensitivity** — refit under an alternative prior; does the conclusion move? (RQ4)

Write diagnostics beside the fit in `outputs/metrics/`. Walkthrough: `docs/bayesian-workflow.md`.
