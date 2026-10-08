# Bayesian workflow for the censored AFT models

`BayesianAFTModel` fits a prespecified Weibull model and supports lognormal family and
wider-prior sensitivity fits. Every fit receives training covariates and outcomes only.
The caller fits imputation, encoding and continuous-variable scaling on train, then passes
the resulting numeric DataFrame. All settings live in `config/default.yaml`.

## Model and assumptions

For Weibull, `T | x, parameters ~ Weibull(shape=k, scale=exp(intercept + x beta))`.
For lognormal, `log(T) | x, parameters ~ Normal(intercept + x beta, sigma)`.
A coefficient exponentiates to a conditional time ratio: values above one are associated
with longer characteristic survival. These are prognostic associations.

**Weibull with a common shape is both an AFT and a proportional hazards model.** Its hazard
ratio for a one-unit predictor difference is `exp(-k beta)`. It does not solve a violation
of proportional hazards merely by using an AFT parameterization. The lognormal model allows
nonproportional hazards but still imposes a shared distributional shape and linear log-time
association. Lognormal is a prespecified sensitivity, never selected using calibration or test.

Right-censored patients contribute the probability of surviving their observed follow-up;
uncensored events contribute a density. `pm.Censored` uses an infinite upper bound for events
and the observed time as upper bound for censored rows. Continuous `pm.Censored` is used only
as the likelihood. Covariates are `pm.Data`; prediction analytically evaluates the stored
posterior and never refits a model. This follows the [PyMC censored distribution API](https://www.pymc.io/projects/docs/en/latest/api/distributions/censored.html).

The likelihood requires censoring to be noninformative conditional on the modeled predictors.
If the endpoint is disease-specific survival, treating other-cause deaths as censored estimates
net survival under an additional competing-event assumption; it is not the observed cumulative
incidence of breast-cancer death. Overall survival has a more direct all-cause interpretation.

## Fixed priors and checks before posterior sampling

The following prior values are fixed protocol choices in months. They are not estimated from
cohort outcomes or from calibration/test data.

| Parameter | Default prior | Wider-prior sensitivity |
| --- | --- | --- |
| Intercept | `Normal(log(120), 0.5)` | `Normal(log(120), 0.75)` |
| Each standardized coefficient | `Normal(0, 0.35)` | `Normal(0, 0.70)` |
| Weibull shape | `LogNormal(log(1.5), 0.35)` | Same |
| Lognormal residual scale | `HalfNormal(0.75)` | Same |

The intercept places the baseline time scale around ten years with substantial uncertainty.
The coefficient prior regularizes log-time associations; doubling its standard deviation is a
meaningfully less regularized sensitivity. These are modeling choices, not estimates of a
population survival distribution. Continuous features are scaled on train; binary indicators
retain their category contrast. Combined predictor variation can still produce long tails,
which makes the prior predictive check essential.

The workflow first samples prior parameter draws and event times at training covariates. It
records time quantiles and the fraction outside configured plausibility bounds. The current
bounds are 0.01 to 1,200 months with at most 20% outside. This is a deliberately broad failure
screen, not a clinical validation of the priors. A failing screen raises before NUTS sampling.
Persisted prior draws and training PPC curves support visual inspection of the resulting time
scale. Do not tune priors against held-out results.

## Sampling and the acceptance gate

All fitting runs in a SLURM allocation. `config.mcmc` supplies chains, tune, draws and target
acceptance; cores are bounded by the allocation. Independent streams derived from the supplied
NumPy Generator control prior parameters, prior event times, NUTS and posterior checks.

Initialization is explicitly `adapt_diag`. Each chain starts at a small, independent,
seeded perturbation around the fixed prior centers; the perturbation is bounded by
`mcmc.initialization_jitter_fraction` (0.1) times the corresponding prior scale. Lognormal
residual-scale initialization uses that fraction on log scale. The model checks both the
log probability and every compiled gradient for finiteness before sampling. These starting
points are numerical settings, not fitted priors. The prior and acceptance gates are unchanged.

This replaces the [PyMC default unit jitter](https://www.pymc.io/projects/docs/en/stable/api/generated/pymc.init_nuts.html)
after a training-only OS molecular fit exposed a numerical failure: one chain never moved,
accepted no proposals, and accumulated 1,000 divergences while its step size collapsed. Its
starting log probability was finite but all 23 compiled gradients were nonfinite; the analytic
Weibull likelihood derivatives remained finite. The censoring tail was near floating-point
underflow. Other chains mixed normally. Prior-centered starts had finite compiled gradients.
The failed fit is retained for diagnosis and excluded from results; all final fits use the
same newly registered initialization protocol.

After NUTS, every stochastic parameter, including every coefficient, must satisfy the configured
gate: maximum rank-normalized R-hat <= 1.01, minimum bulk and tail ESS >= 400, and zero
divergences. Nonfinite diagnostics also fail. `BayesianDiagnosticsError` prevents prediction
and downstream reporting from a failed fit. The caller may save the failed inference object
and diagnostics for investigation before stopping the stage.

Inspect convergence or geometry failures before changing settings. More draws address Monte
Carlo precision when chains otherwise mix; they do not fix nonidentifiability or a poorly
specified model. Raising target acceptance alone does not establish good geometry.

PyMC 6 and ArviZ 1 use `xarray.DataTree` for inference storage, replacing the older
`InferenceData` class. `.nc` artifacts retain posterior, sample statistics, prior and data
groups, while a JSON sidecar retains config, feature names, diagnostics and aggregate workflow
summaries. Load calls recheck diagnostics. See [ArviZ's storage schema](https://python.arviz.org/en/stable/schema/schema.html).

## Training posterior predictive check

The check reports the training Kaplan-Meier estimate and its log-log confidence interval on
the configured grid, strictly below maximum training follow-up. Posterior survival curves
average over the empirical training predictor distribution, with a pointwise posterior band.
The corresponding prior survival band shows how learning changed uncertainty. These bands
measure uncertainty about the survival curve; they are not individual event-time intervals.

A second check generates latent event times and independently generates censoring times from
the training reverse Kaplan-Meier distribution. It compares the replicated event fraction to
the observed training event fraction. Unidentified censoring-tail mass stays at infinity. This
is a descriptive check under marginal independent censoring, not proof of conditional
independence and not evidence that the outcome model learned the censoring mechanism.

`workflow_["training_ppc"]` stores these curves and absolute discrepancies. Model inadequacy
must be reported even when the sampler passes. The fixed Weibull/lognormal comparison and
wider-prior refits are orchestrated by stage 03 for each feature set using identical training
rows. Compare training PPC, coefficient summaries and posterior prediction changes; never
choose a sensitivity variant based on calibration/test performance.

## RQ4 and prediction

`predict_survival_bands` returns pointwise parameter-uncertainty bands. `predict_quantiles`
inverts the posterior mixture of event-time distributions, combining parameter uncertainty
with future event-time variation. It does not average component quantiles: that would be a
different distribution. A fixed evenly spaced subset of posterior draws controls prediction
cost and is recorded by `bayesian.predictive_draws`.
Coefficient means and credible intervals use all retained posterior draws.

Mixture quantiles are solved on log time between the minimum and maximum component quantile,
without an arbitrary time cap. The reported predictive interval and a conformal interval
answer different questions. The former depends on priors and the parametric model; empirical
conformal coverage and width are evaluated separately using the frozen partitions.

Scientific results belong in `outputs/metrics/` and are written only after actual fits,
independent review and recorded manifests. Synthetic smoke-test diagnostics verify the
implementation and are never METABRIC results.
