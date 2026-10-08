# Registered analysis decisions — 2026-09-29

The user selected breast-cancer-specific survival (DSS) as primary and overall survival
(OS) as sensitivity, and delegated the remaining choices. These decisions precede
calibration and test inspection. Both endpoints and feature arms use the same DSS-event
stratified split: 989 training, 495 calibration and 495 test patients.

The public cBioPortal clinical extract yields 1,979 eligible patients after excluding
528 without survival follow-up, one without vital-cause classification, and one with
nonpositive follow-up. These are sequential exclusions. Numeric predictors retain
missingness: medians, scales and categories are fitted on training patients only.
Missing indicators accompany numeric imputation; unknown categorical levels have an
explicit encoding. This retains 165 patients who would be excluded by complete-case
analysis. It is a reproducible prediction strategy, not a claim that missingness is random
or a substitute for multiple-imputation inference about coefficients. See the
[prediction-time missingness study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7586995/).

Retain the supplied claudin-low label separately and map NC/missing subtype to unknown.
Claudin-low is a recorded phenotype, not asserted to be a sixth intrinsic PAM50 subtype;
its relationship to intrinsic subtypes is discussed by
[Fougner et al.](https://www.nature.com/articles/s41467-020-15574-5).

Retain full recorded follow-up. Thirteen eligible patients have follow-up beyond
300 months, including two DSS events and five OS events. A blanket truncation would
discard information. Late recurrence is relevant to METABRIC, as documented by
[Rueda et al.](https://www.nature.com/articles/s41586-019-1007-8).
Prediction-metric evaluation uses a training-derived interior grid capped at 120 months;
this does not alter fitting labels. Conformal prediction targets the restricted time
min(T, 300 months), with a prespecified 120-month sensitivity analysis. Restricted-time
coverage is not coverage of the unrestricted event time.

DSS treats other-cause death as censoring. Its survival function describes net
cause-specific survival under the model's censoring assumptions, not crude cumulative
breast-cancer mortality in the presence of competing deaths. Marginal reverse-KM IPCW
also assumes censoring independence from the event time and relevant predictors;
violations can affect performance metrics and IPCW coverage estimates. Observable
coverage bounds and the conservative conformal construction are reported separately.

Primary models are Cox PH, Bayesian Weibull AFT and untuned RSF. Bayesian lognormal and
wider-prior Weibull fits are training-only sensitivity analyses, not selected using
calibration or test performance. Primary ranking uses Uno C, with IBS and uncertainty
width/coverage presented alongside it. Bootstrap intervals condition on the frozen fits
and calibration; they do not include model-fitting or split-selection uncertainty.

Every compute run uses SLURM. Explicitly allowed uncommitted runs preserve a deterministic,
content-addressed source/config archive alongside the manifest; Git SHA, dirty state,
config/artifact hashes, dependency versions and job ID are retained. Source snapshots
make these runs traceable before the user creates checkpoint commits. Reproduction from
a commit alone requires that commit's source to match the recorded archive.

Training diagnostics detected global proportional-hazards departures in all four Cox fits.
They are retained as registered prognostic benchmarks, with violations reported. Weibull
AFT also implies PH; the lognormal training sensitivity and PPCs assess this limitation
without choosing a family using held-out outcomes. Unrestricted net-DSS predictive tails
can exceed plausible lifetimes: they are parametric extrapolations, not literal clinical
forecasts. The registered restricted-time conformal targets avoid interpreting those tails
as empirically supported unrestricted event times.
