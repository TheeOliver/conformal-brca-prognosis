# Classical survival models

`CoxPHModel` and `RandomSurvivalForestModel` implement the shared model interface on
a numeric pandas DataFrame. Preprocessing must be fitted exclusively on training
data before calling either adapter. Both feature sets use the same frozen patient
split. Every model validates the original feature names and their order at prediction.
The outcome is a structured array with a boolean event indicator and positive,
finite survival time in months. Censored patients remain in the fit.

## Cox proportional hazards

`CoxPHModel(config.cox, rng=rng)` uses scikit-survival's Cox estimator. The config
supplies `alpha`, `ties`, `n_iter`, `tol`, `ph_diagnostic_transform` and
`ph_diagnostic_alpha`. Its default unpenalized fit uses Breslow ties; Efron ties
are supported too. Optimization warnings terminate fitting rather than publishing
unconverged predictions. Coefficients describe conditional prognostic associations.
[Estimator documentation](https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.linear_model.CoxPHSurvivalAnalysis.html)

Features that are constant in the training matrix are omitted from the Cox fit and
listed in `constant_features_`. This handles explicit unknown-category columns
that contain no training observations. Their association is unidentified and the
adapter assigns them no learned contribution at prediction. Nonconstant collinear
columns are rejected: the preprocessor must omit one reference category per factor.
`active_feature_names_` maps fitted coefficients back to retained feature names.

The adapter computes `ph_diagnostic_` from the training data alone. It tests adding
time interactions to the fitted linear predictor. With configured transformation
`g(t)` (log time by default), the alternative is `beta(t) = beta + gamma*g(t)`.
The score and information come from the same Breslow or Efron partial likelihood
as the fit; the covariance accounts for estimation of the ordinary Cox coefficients.
The information blocks and score construction follow the maintainer's description
of the score test in the [survival vignette, Computational details](https://therneau.r-universe.dev/survival/doc/survival.pdf).

The stored result contains one global chi-square statistic and one per retained
encoded column, with degrees of freedom, p-values, event count and the configured
flagging level. Per-column p-values are unadjusted exploratory diagnostics. A
categorical predictor is represented by its encoded-column tests rather than a
pooled factor test. The reference distribution is asymptotic. Failure to reject
does not establish proportional hazards; rejection records a model limitation and
does not trigger refitting or selecting a model on calibration/test data. A
penalized fit or singular information yields an explicit unavailable diagnostic,
never invented p-values. The approach agrees with the role of score tests described
in the [survival cox.zph documentation](https://therneau.r-universe.dev/survival/doc/manual.html#cox.zph).

Only aggregate diagnostic results are retained. The tests compare the statistic
against independent numerical derivatives of the time-interaction likelihood,
including tied event times, and check detection of a known synthetic violation.

## Untuned random survival forest

`RandomSurvivalForestModel(config.rsf, rng=rng)` supplies the configured
`n_estimators`, `min_samples_split`, `min_samples_leaf`, `max_features`, and `n_jobs`
to scikit-survival. There is no tuning. A seed is drawn once from the injected
generator and retained, so repeated fits of the adapter reproduce the same forest.
`n_jobs` is a positive CPU count within the SLURM allocation; `-1` is rejected.
Risk predictions have the estimator's convention that larger values mean worse
prognosis. Survival curves average the trees' terminal-node Kaplan–Meier curves.
[Estimator documentation](https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.ensemble.RandomSurvivalForest.html)

## Survival steps, tails and quantiles

Both adapters return survival probabilities on the caller's grid as right-continuous
step functions. Survival is one before the first observed time; at a knot the value
after that jump is used. Beyond the final training knot the last probability is
carried forward. That computational convention adds no information about the
unobserved tail. Evaluation grids must remain within training follow-up support.

`predict_quantiles(X, probabilities)` returns event-time quantiles for probabilities
strictly between zero and one. A quantile is the first fitted knot where the
estimated event CDF reaches the requested probability. It is positive infinity
when the curve never reaches that probability. In particular, the final observed
follow-up time must never be substituted for an unidentified median or upper bound.
Consumers must retain or explicitly report these infinite values rather than
silently truncating them.

`save(path)` and the class method `load(path)` persist the estimator, feature
contract, random seed/state and Cox diagnostic in a versioned pickle artifact.
Only load trusted project artifacts. They belong in gitignored `outputs/models/`;
they are generated fits and must not enter git.
