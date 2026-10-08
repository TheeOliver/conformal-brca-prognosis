# Conformal prediction for survival data

> Implementation: `src/brca/conformal/survival.py`. All scores and rank thresholds are
> hand-written; reverse Kaplan–Meier estimation uses scikit-survival.

The post hoc [uncertainty reliability experiment](uncertainty-reliability.md) adds
conservative one-sided lower bounds across the three model families and a fixed-ridge
Qin sensitivity with a strict all-refits-successful gate. Its outputs are isolated from
the primary analysis below. Lower bounds answer a different prediction question; ridge
addresses numerical stability without repairing censoring support.

## Estimand and declared analyses

The primary interval predicts **restricted individual survival time**
`Z = min(T, tau)`, in months, with the prespecified horizon `tau = 300` months.
The prespecified secondary horizon is **120 months**, processed by the same
functions and declared in `conformal.horizons_months: [300, 120]`, with
`conformal.primary_horizon_months: 300`. Both are reported; neither is selected
from held-out outcomes, censoring support or achieved coverage.
An interval ending at 300 includes patients surviving beyond 300; it does **not**
assert that their event occurs by 300. This restriction applies to the conformal
estimand independently of any administrative censoring used for model fitting.
Predictions use `alpha = 0.20, 0.10, 0.05` for 80%, 90%, 95% nominal coverage.

Disease-specific survival is the primary endpoint, with overall survival as a
sensitivity on the **same patient split**. For DSS, censoring deaths from other
causes defines a net survival analysis: a latent time to breast-cancer death
compatible with survival until the competing death. That latent time is not
observed after a competing death. It does not identify crude breast-cancer
mortality or its cumulative incidence. Interpreting the fitted DSS survival
model or IPCW estimate as net survival requires assumptions about independence
of competing deaths and breast-cancer survival. The conservative interval
construction below only uses the information `T > censoring_time` under this
latent-time interpretation; it cannot validate that interpretation.

The primary method is **conservative score-envelope split conformal prediction**.
An optional, prespecified **IPCW score sensitivity** explores efficiency under
stronger assumptions. The two methods are reported together without choosing the
winner after looking at test coverage. Neither subgroup analyses nor different
alpha levels are used to select or refit a score.

## Primary score and calibration

Let `[l(x), u(x)]` be a base model's predictive quantile interval, fitted only on
train and clipped to `[0, tau]`. The same nonnegative conformalized quantile
regression score is used for every model:

```text
s(x, z) = max(l(x) - z, z - u(x), 0).
```

We observe `Y = min(T, C)` and `event = 1{T <= C}`. If the event is observed or
`Y >= tau`, the restricted target `Z = min(Y, tau)` is known. Otherwise its
compatible values are `(Y, tau]`. Use the following calibration score envelope:

```text
known restricted target:  s_bar = max(l - Z, Z - u, 0)
censored before tau:      s_bar = max(l - Y, tau - u, 0).
```

All censored patients remain in calibration. Their censoring times are never
silently treated as event times. The score is nonnegative, so calibration can
widen a base interval but cannot shrink it. This sacrifices some efficiency
when a base model already gives intervals wider than needed.

Within a calibration sample of size `n`, compute exactly

```text
k = ceil((n + 1) * (1 - alpha))
q = sorted(s_bar)[k - 1] if k <= n else +infinity
C(x) = [max(0, l(x) - q), min(tau, u(x) + q)].
```

There is no interpolation and no clipping of the rank to `n`. A missing rank
means an unbounded pre-intersection set. For restricted survival its intersection
with the known target support is `[0, tau]`: the report records both an infinite
calibration threshold and a full-support interval. Finite restricted width must
not conceal this loss of information. For unrestricted survival the conservative
fallback assigns infinity to censored scores, often producing `[0, infinity)`;
its infinite width is explicitly reported.

## Why the conservative construction covers

Conditional on the training fit, suppose calibration and a new patient are
exchangeable. Their observed score envelopes are exchangeable. The ordinary
split-conformal rank argument gives

```text
P(s_bar_new <= q) >= 1 - alpha.
```

Every latent target compatible with an observation has
`s(x, Z) <= s_bar`. Therefore the event `s_bar_new <= q` implies
`Z_new in C(X_new)`, giving the same lower coverage bound. Independent censoring
is **not** required for this domination argument. It requires the observation
interval to contain the intended latent outcome and the calibration/new-patient
exchangeability assumption to hold. It promises marginal coverage over the
calibration and future-patient draw; it does not promise an exact coverage level,
coverage conditional on a fixed realised calibration set, or coverage for each
PAM50 or ER subgroup.

This is a conservative application of interval-outcome conformal inference.
Liu, de Paula and Tamer develop prediction sets under interval censoring and an
interval quantile-regression construction in Appendix E. This implementation
uses the simple nonnegative score above and does not claim their efficiency
results or reproduce their optimal-set estimator.
[Primary paper](https://arxiv.org/abs/2501.10117).

### The frozen event-stratified split

The project's split fixes event counts in each partition. Consequently, blindly
applying an iid marginal rank theorem to all calibration rows is not justified.
The primary method calibrates separately in each **primary DSS event stratum**,
then uses the **maximum** stratum threshold for every prediction. This procedure
does not need a future patient's event label. Within-stratum exchangeability
gives the preceding rank argument in each stratum; taking the maximum contains
every stratum-specific interval, hence also covers any mixture of those strata.
A required stratum absent from calibration contributes an infinite threshold.

The OS sensitivity must still pass the **primary DSS event stratum** to
calibration, because that variable determined the shared split. Do not replace
it with the secondary endpoint's event indicator. The metadata records each
stratum's calibration count and threshold, plus total `n_cal`.

Stratum exchangeability remains an assumption: secular treatment changes,
non-random sampling, or preprocessing applied differently between partitions can
break it. The code does not establish that real METABRIC patients satisfy it.

## Optional approximate IPCW score sensitivity

For this sensitivity, estimate the censoring distribution on **train only** using
reverse Kaplan–Meier with event-before-censor handling of tied times. Write
`G_ge(t) = P(C >= t) = G(t-)`. Restricted targets observed as an event before the
horizon, or known to survive to the horizon, have weights `1/G_ge(Z)`; targets
censored earlier contribute no observed target score and have weight zero. All
rows still contribute to cohort counts and diagnostics.

Take the `(1-alpha)` weighted quantile of the observed target scores and an
infinity atom of weight `1/G_ge(tau)`. With no censoring and unit weights this
reduces exactly to `ceil((n+1)(1-alpha))`. Estimated outcome-dependent weights
are not equivalent to ordinary unweighted exchangeable scores. This is a
**sensitivity method without an asserted finite-sample guarantee**. It is pooled;
the project's event-stratified sampling creates an additional reason not to
assert its ordinary iid weighted-rank guarantee.

Marginal censoring weights require `C` independent of **both** `(T, X)`, not merely
`T independent of C given X`. A marginal Kaplan–Meier curve cannot correct
covariate-dependent censoring by itself. In DSS this assumption includes other
causes of death and is especially strong. Estimates require positive support
through `tau`; if train follow-up does not reach `tau` or `G_ge(tau)` is below the
configured minimum, calibration returns the whole target support and records why.
Weights are never silently clipped and unsupported tails are never extrapolated.

Candès, Lei and Ren study lower predictive bounds under Type-I censoring, with
conditions specific to that observation scheme. Their finite-sample result
cannot be transferred directly to two-sided intervals in general right-censored
METABRIC data. [Primary paper](https://doi.org/10.1093/jrsssb/qkac004).
Qin and colleagues propose bootstrap conformal intervals for general censoring;
their results highlight the censoring-model and upper-tail difficulties.
The primary registered score envelope here is distinct from Qin's method. A
post hoc, two-sided Cox working-model implementation of their bootstrap CPI is
documented in [the Qin exploratory report](qin-exploratory.md); its results are
kept separate from the frozen primary analysis. See also the local paper
`references/literature/qin_2024_conformal_predictive_intervals_survival.md`.

## What can be measured on censored test data

Latent coverage is not directly observed for patients censored before the
horizon. The report first gives assumption-light **observable coverage bounds**:

- A known target is certainly/possibly covered exactly when it lies in its interval.
- A censored target is certainly covered when the entire compatible set
  `(Y, tau]` lies inside the interval: `lower <= Y` and `upper >= tau`.
- It is possibly covered when the interval intersects that set:
  `upper > Y` and `lower <= tau`.

The mean of the first indicator lower-bounds realised latent test coverage; the
mean of the second upper-bounds it. Wilson intervals describe sampling
uncertainty in these **observable indicator rates**, not a binomial confidence
interval for the unknown latent coverage. The same bounds are calculated for
Bayesian predictive intervals and other uncalibrated model intervals.

The assumption-dependent IPCW coverage estimate is a self-normalised (Hájek)
weighted fraction of covered known targets. A Horvitz–Thompson version divides
the weighted covered sum by total test `n` and is supplied as a diagnostic; it
can exceed one in a finite sample. The report also supplies the estimated total
observation mass, maximum weight and effective sample size. Agreement of these
estimators is not evidence that the censoring assumptions are true.

Percentile bootstrap confidence intervals resample **test patients**, keeping
model predictions, calibration and the train-fitted censoring distribution
fixed. They describe test-sample uncertainty conditional on those fits and omit
training/calibration/nuisance-estimation uncertainty. The same seed and aligned
patient order permit paired resampling across model/feature arms. No binomial
confidence interval is attached to weighted IPCW coverage. Unsupported censoring
or empty subgroups produce explicit null estimates and status labels, not zeros.

Pass `bootstrap_strata` as the **primary DSS event indicator** for the original
ordered test rows. Resampling takes place within these strata, keeping their
sample counts fixed as in the frozen split. For subgroup reports the same strata
are restricted to the subgroup first. This applies to both DSS and OS; the OS
event indicator must not replace the variable that originally defined the split.
Ordinary patient bootstrap remains available when no strata are supplied, for
settings without a stratified design. The resampling method is recorded.

Each bootstrap statistic requires at least **two** defined replicates and the
configured `minimum_bootstrap_success_fraction` (protocol value 0.90). Success and
failure counts accompany every CI; samples are not redrawn. Any CI retained after
discarding undefined replicates is explicitly labelled
`conditional_on_successful_resamples` and approximate. The success-fraction gate
does not restore unconditional confidence coverage. Below the gate the interval
is null with `bootstrap_gate_failed`, while its point estimate remains visible.

The same bootstrap draws also provide intervals for both observable coverage
bounds, the median and quartiles of width, full-support fraction and infinite
interval fraction. These observable statistics retain their uncertainty
intervals even when IPCW estimation is unsupported. Width summaries and their
bootstrap limits use empirical order-statistic quantiles to handle infinite
widths without invalid interpolation. Wilson intervals are additional descriptive
approximations for the observable rates, not exact guarantees for a stratified
test sample.

For the paired feature-set comparison, `paired_interval_comparison` uses the
**identical ordered test patients and identical bootstrap selections** for the
clinical and clinical+PAM50 arms. It reports molecular-minus-clinical differences
in IPCW coverage, both observable coverage bounds, median/quartiles of width,
full-support and infinite-interval fractions. A positive width difference means
wider molecular intervals; the statistic subtracts the two marginal medians,
not the median of each patient's width difference. Models and calibration remain
fixed during this paired test bootstrap. An undefined infinity-minus-infinity
width difference is labelled unavailable, never replaced with zero.

## Efficiency and subgroup reports

Every coverage report includes calibration and test counts, median and quartiles
of width, proportion of infinite intervals and proportion covering the full
known support. An all-support interval can achieve perfect coverage while being
uninformative. This is a substantive result, especially under heavy DSS censoring.
Infinite widths are encoded as JSON null with explicit infinity flags, never as
non-standard JSON `Infinity` or `NaN` values.

PAM50 and ER reports reuse the marginal intervals and show subgroup sample sizes.
They measure variation in coverage; they are not conditional-coverage guarantees.
Small groups yield unstable IPCW estimates and broad uncertainty; missing levels
are retained as `n_test=0` with null estimates. Subgroup findings should be
presented with these limitations and without causal language.

## Validation contract

Synthetic tests cover the corrected finite-sample rank and too-small-calibration
case; known and censored score envelopes; repeated independent calibration/test
simulations with uncensored near-nominal coverage; conservative latent coverage
under deliberately dependent censoring and event-stratified allocation; heavy
censoring, full-support and infinite intervals; censoring-time ties and support;
IPCW estimation/calibration under independent censoring; restricted horizon
observability; missing subgroups; and strict JSON serialization.

Synthetic validation establishes implementation behaviour under the simulated
conditions. It is not evidence of achieved coverage on METABRIC. A separate
conformal-validity review must examine the implementation before real coverage
numbers are presented.
