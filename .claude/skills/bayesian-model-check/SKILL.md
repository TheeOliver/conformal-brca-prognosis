---
name: bayesian-model-check
description: Run the full Bayesian workflow on a PyMC survival fit - prior predictive check, sampling, R-hat/ESS/divergence gate, posterior predictive check and prior sensitivity. Use when asked to check a Bayesian fit, investigate divergences or bad R-hat, justify priors, run posterior predictive checks, test prior sensitivity, or report what uncertainty the Bayesian model gives (RQ4).
---

# Check a Bayesian fit

This answers RQ4 — *what uncertainty does the probabilistic model itself give?* — and it is
also the gate that decides whether a fit may be reported at all.

Read `docs/bayesian-workflow.md` before starting. Thresholds live in
`config.mcmc.diagnostics`, never hard-coded.

## 1. Prior predictive check — before any sampling

Draw from the prior predictive and plot simulated survival times.

Ask: are these times *clinically possible*? An AFT log-scale prior that looks harmlessly wide
routinely implies survival of 10^5 months. If the prior predictive is absurd, fix the prior
now — fitting first and discovering this later wastes the whole run.

## 2. Sample

Settings come from `config.mcmc`. Sampling is compute — submit via `sbatch`, never on the
login node. Persist the result as InferenceData to `outputs/models/<name>.nc`.

## 3. The diagnostics gate

Assert, and **raise** on failure:

| Check | Threshold |
| --- | --- |
| max R-hat | <= 1.01 |
| bulk ESS | >= 400 |
| tail ESS | >= 400 |
| divergences | 0 |

If it fails, diagnose rather than resampling blindly:

- **Divergences** → inspect the pair plot for the funnel; try a non-centred parameterisation,
  then raise `target_accept`. Raising `target_accept` alone masks the geometry problem.
- **Low tail ESS, fine bulk ESS** → the posterior mean is fine but the interval endpoints are
  noise. Since this project reports intervals, this is a blocking failure, not a nuisance.
- **R-hat > 1.01** → check for an unidentified parameter or a label-switching symmetry before
  simply sampling longer.

## 4. Posterior predictive check

Overlay posterior predictive survival curves on the **training** Kaplan-Meier, with the KM
confidence band. Systematic deviation in the tail is the one that matters here: it is exactly
the region where individual prognosis is hardest and where conformal intervals get widest.

Also check the posterior predictive censoring pattern is plausible — a model that reproduces
event times but not the censoring structure has learned the wrong thing.

## 5. Prior sensitivity

Refit under at least one genuinely different prior (e.g. tighter scale, heavier tails) and
report whether conclusions move. With ~2000 patients the likelihood should dominate; if it
does not, that is a finding about how much of the reported uncertainty is prior, not data —
and it belongs in the Discussion.

## 6. Record

Write every diagnostic to `outputs/metrics/bayes_diagnostics_<name>.json` so a reviewer never
re-runs sampling to check. Report to the user: gate pass/fail per check, what the PPC showed,
and whether priors mattered.
