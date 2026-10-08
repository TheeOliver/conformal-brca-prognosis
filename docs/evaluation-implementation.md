# Evaluation implementation

`src/brca/evaluation/metrics.py` evaluates fixed model predictions against survival
labels supplied by stage 05. It loads no files and never refits a predictive model.
Only stage 05 reads the frozen test split. Inputs are structured boolean-event /
floating-point-months arrays and a mapping from model identifiers (for example,
`cox__clinical` and `cox__clinical_molecular`) to `risk` and `survival` arrays.

## Shared times and censoring support

Call `training_time_grid(y_train, config)` before predicting. It computes configured
percentiles of observed training follow-up, caps them at the configured maximum
horizon and strictly before the last training event, removes duplicate knots and
requires at least two times strictly above the minimum training follow-up. The
default horizon is 120 months. `evaluate_models` verifies that its supplied grid
equals this training-only construction. Test outcomes cannot move the grid.

The full held-out sample must contain an observed event by the first grid time,
have follow-up below that time, and retain observations beyond the final grid time.
Failure raises `EvaluationSupportError` carrying aggregate diagnostics. Training
reverse Kaplan–Meier censoring survival must be positive everywhere used for
evaluation. Undefined full-sample metrics also stop evaluation rather than
silently dropping a model or choosing new horizons.

For evaluation only, follow-up after the final grid time is administratively
censored at the next representable float above that time. Events at the final
grid time remain events; later observations remain controls at every grid time.
This avoids asking the training censoring estimator to extrapolate for irrelevant
late follow-up. The returned report records the exact cap and number of affected
observations. The input labels are copied, so fitted model targets and original
follow-up remain intact. The same cap is used in every bootstrap replicate.

## Metrics and interpretation

Uno's concordance uses the scalar `risk` prediction, with larger risk indicating
worse prognosis, and a truncation time immediately above the last grid knot.
Censoring weights are always estimated from the training labels, following
[scikit-survival's IPCW concordance API](https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.concordance_index_ipcw.html).

The time-specific Brier score and integrated Brier score use survival probabilities
on the same shared grid. Lower scores indicate more accurate survival predictions;
they combine calibration and discrimination. The implementation calls the public
[Brier score](https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.brier_score.html)
and [integrated Brier score](https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.integrated_brier_score.html)
functions directly.

Cumulative/dynamic AUC uses the time-dependent event risk `1 - S(t)` at each grid
knot. Thus it measures ranking at each horizon, including models whose ordering
can change over time. Both the AUC curve and the survival-weighted mean AUC returned
by the library are reported. Higher values indicate better discrimination.
[AUC API](https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.cumulative_dynamic_auc.html)

These marginal IPCW implementations assume censoring is independent of both the
event time and predictors: `C independent of (T, X)`. A marginal censoring curve
does not address censoring that depends on included or omitted predictors.
Including age or molecular subtype in the survival model does not make a marginal
censoring estimator adjust for those variables. Conditional independence given
predictors alone is insufficient for the weights used here. That limitation
accompanies the reported performance measures, especially for net DSS analyses
that censor other-cause deaths.

## Paired uncertainty

The configured number of bootstrap samples is drawn over test patients. Stage 05
passes `bootstrap_strata` containing the aligned boolean **primary DSS event**
indicators that governed the frozen split, including when evaluating OS. Patients
are resampled with replacement within each original stratum, preserving its count.
One resulting index vector per replicate is applied to every model and feature set. Every
scalar and every curve point receives a percentile interval at the configured
confidence level. Training censoring estimates and fitted predictions are held
fixed; these intervals describe held-out sampling uncertainty conditional on the
training fit and fixed stratum counts, and do not include variation from refitting
on new training cohorts. Omitting `bootstrap_strata` requests an ordinary bootstrap;
the output identifies it as an IID approximation for the stratified split because
it additionally varies the stratum composition.

Feature comparisons subtract clinical from clinical plus molecular metrics within
each replicate. They report paired intervals for Uno C, IBS, mean AUC and every
Brier/AUC knot. Positive differences favor the molecular model for Uno/AUC; negative
differences favor it for Brier/IBS. Identical predictions must produce exact zero
paired differences and intervals; this is checked in the tests, along with a
manual bootstrap comparison using nonidentical survival curves.

A bootstrap sample may have no cases, no comparable pairs, or insufficient follow-up
for a requested metric. Such metric estimates are explicitly counted as failed,
with no redraw. Each summary records successes, failures and the success fraction.
If fewer than the configured fraction (default 90%) succeed, or fewer than two
estimates exist, that interval is null and its status is `bootstrap_gate_failed`.
The whole report then has `reportable: false`. When the gate passes but any draws
were undefined, the interval, model and overall report have status
`conditional_approximate`: the interval conditions on a defined estimate rather
than representing the complete bootstrap distribution. Such an interval remains
reportable only with this qualification and its failed-draw count. It is never
labelled unqualified `ok`. No NaN or infinity is serialized.
Bootstrap intervals are exploratory and have no multiplicity adjustment; the
configured primary metric is recorded to avoid selecting a favorable metric later.

## Output contract

`evaluate_models(y_train, y_test, predictions, times, config, rng,
bootstrap_strata=primary_dss_event)` returns a JSON-safe
dictionary with `models`, `paired_feature_comparisons`, shared-grid and censoring
metadata, and bootstrap settings. Each model contains `uno_c`, `ibs`, `mean_auc`,
`brier` and `auc`. Scalars have `estimate`, `ci_lower`, `ci_upper`, bootstrap counts
and status; curve entries additionally have `time_months`. No patient rows or
per-patient predictions are included in this report.
