# Bayesian workflow for the Weibull / AFT model

> Pulled in by: `bayesian-model-check` skill; `.claude/rules/bayesian-pymc.md`,
> `.claude/rules/mcmc-diagnostics.md`.

## Model

Accelerated failure time: covariates scale time directly, so `exp(beta)` is a **time ratio**,
easier to explain to a clinician than a hazard ratio — and, unlike Cox, it needs no
proportional-hazards assumption.

```
log T_i = mu + x_i . beta + sigma * eps_i
```

Weibull baseline. Right censoring enters through `pm.Censored`, the single censoring idiom in
this project.

## Priors, and why they are not innocuous

Standardise continuous predictors on train first, so one prior scale means the same thing for
every covariate.

| Parameter | Prior | Reasoning |
| --- | --- | --- |
| `mu` (intercept) | `Normal(log(median_followup), 1)` | centred on the observed time scale |
| `beta` | `Normal(0, 0.5)` | on standardised covariates: time ratios mostly within ~2.7x |
| `sigma` / shape | `HalfNormal(1)` | positive, weakly informative |

**A prior that looks vague on the log scale is not vague on the time scale.** `Normal(0, 10)`
on `mu` puts real mass on survival times of `exp(10)` months. This is the most common way to
get an absurd prior predictive here — which is why the prior predictive check comes first.

## Order of operations

1. **Prior predictive check.** Draw and plot simulated survival times. Clinically possible?
   Fix the prior now, before spending a sampling run.
2. **Sample.** `config.mcmc`: 4 chains, 1000 tune, 1000 draws, `target_accept` and
   `random_seed` from config. Submit via `sbatch`.
3. **Diagnostics gate.** R-hat <= 1.01, bulk and tail ESS >= 400, zero divergences. Assert and
   raise. A failing fit is not a result.
4. **Posterior predictive check.** Overlay posterior predictive survival on the training
   Kaplan-Meier with its confidence band. Watch the tail: that is where individual prognosis is
   hardest and where conformal intervals widen most.
5. **Prior sensitivity.** Refit under a genuinely different prior. With ~2000 patients the
   likelihood should dominate; if it does not, report how much of the stated uncertainty is
   prior rather than data.

## Diagnosing failures

| Symptom | Likely cause | First move |
| --- | --- | --- |
| divergences | funnel geometry | non-centred parameterisation, *then* raise `target_accept` |
| low tail ESS, fine bulk ESS | interval endpoints unresolved | blocking here, since this project reports intervals |
| R-hat > 1.01 | unidentified parameter or multimodality | check identifiability before sampling longer |
| very slow sampling | unstandardised covariates | standardise on train |

Raising `target_accept` to silence divergences without fixing the geometry hides the problem
rather than solving it.

## What this answers

RQ4 — the uncertainty the probabilistic model reports *given its priors and structure*. That is
a different claim from conformal coverage, which is checked empirically against held-out data.
Keeping the two apart is central to the thesis.

## TODO

- [ ] Record the final priors and the justification for each.
- [ ] Record prior-sensitivity findings.
