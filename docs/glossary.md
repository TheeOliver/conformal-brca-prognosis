# Glossary

> Referenced from CLAUDE.md. The vocabulary this thesis spans three fields, and the same word
> sometimes means different things in each.

## Survival analysis

**Survival function** `S(t) = P(T > t)` — probability of still being event-free at `t`.

**Hazard function** `h(t)` — instantaneous event rate at `t` given survival to `t`. Not a
probability; it can exceed 1.

**Right censoring** — the event had not occurred when observation ended. The patient
contributes "event-free up to `t`", which is information, not a missing value. Dropping
censored patients biases survival downward.

**Administrative censoring** — censoring caused by the study's end date rather than by the
patient leaving. Usually independent of prognosis, which is why the horizon must be recorded.

**Competing risk** — a different event that makes the event of interest impossible (death from
another cause). Treating it as ordinary censoring assumes it is uninformative, which is often
false.

**Kaplan-Meier** — non-parametric estimator of `S(t)` handling censoring; the reference curve
every model should be able to reproduce on its training data.

**Hazard ratio** — ratio of hazards between two covariate values, assumed constant over time
under Cox. Here it is a **conditional association, not an effect**.

**Proportional hazards** — Cox's assumption that covariates scale the hazard by a
time-constant factor. Check it; report violations.

**AFT (accelerated failure time)** — covariates scale *time* rather than hazard. `exp(beta)`
is a time ratio. No proportional-hazards assumption.

## Bayesian inference

**Prior** — the distribution over parameters before seeing data. Never innocuous: a prior that
looks vague on a log scale can be absurd on the time scale.

**Posterior** — parameter distribution after conditioning on the data.

**Posterior predictive distribution** — distribution over *new observations*. This, not the
posterior over parameters, is what gives an individual patient's predictive interval.

**Credible interval** — interval containing a parameter with stated posterior probability.
Distinct from a confidence interval, and distinct again from a conformal interval.

**R-hat** — between- vs within-chain variance. `> 1.01` means chains have not mixed.

**ESS (bulk / tail)** — effective sample size. Bulk governs the posterior mean; **tail governs
interval endpoints**, which is what matters here.

**Divergence** — a failed NUTS trajectory, signalling unexplored geometry and a biased
posterior. Zero tolerated.

**Posterior predictive check** — simulate data from the fitted model and compare to observed.

## Conformal prediction and trustworthy ML

**Conformal prediction** — wraps any model to produce sets with a finite-sample coverage
guarantee, given exchangeability.

**Split / inductive CP** — fit on train, score on calibration, predict on test. The only
variant used here.

**Nonconformity score** — how unusual an observation is under the model. Larger = worse fit.

**Calibration set** — held-out data used *only* to convert scores into a threshold. Spent the
moment it influences any modelling choice.

**Exchangeability** — joint distribution invariant to reordering. Weaker than i.i.d., and the
condition the guarantee rests on. Censoring strains it.

**Marginal coverage** — `P(Y in C(X)) >= 1 - alpha` averaged over patients. What split CP
promises.

**Conditional coverage** — coverage holding *within* a subgroup. **Not** promised; measuring
where it fails is a result.

**alpha** — miscoverage. `alpha = 0.10` gives a 90% interval.

**Validity** — does empirical coverage match nominal.

**Efficiency** — how narrow the interval is. A valid interval spanning all follow-up is useless.

**Calibration (ML sense)** — do predicted probabilities match observed frequencies. Distinct
from *calibration set* above; both words appear in this project.

**Discrimination** — can the model rank patients correctly. Independent of calibration.

**Distribution shift** — test data drawn differently from training data; breaks exchangeability
and so breaks the guarantee.

## Bioinformatics

**PAM50** — 50-gene expression classifier assigning a breast tumour to an intrinsic molecular
subtype: LumA, LumB, Her2-enriched, Basal-like, Normal-like.

**ER / PR / HER2 status** — receptor status from immunohistochemistry; standard clinical
predictors, partly overlapping with what PAM50 captures.

**Batch effect** — technical variation from processing rather than biology; can masquerade as
signal.

**Feature-selection leakage** — choosing features using the whole dataset, then "validating" on
part of it. Classic in gene-expression work, and a reason selection happens inside train only.

## Study design

**Prognostic model** — predicts an outcome. Makes no claim about what would happen under an
intervention. **This project is prognostic.**

**Causal effect** — what would change under an intervention. Not identified by this design.
A strong predictor need not be a cause.
