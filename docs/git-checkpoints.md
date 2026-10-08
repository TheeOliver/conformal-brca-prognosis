# Git checkpoints for the completed thesis analysis

These commands are drafts for the user; no commits, pushes or PR creation were performed.
All checkpoints use the existing `feat/thesis-pipeline` branch and can form one cumulative
PR. Run uncommitted checkpoints in order. Generated data, posterior draws, caches, figures
and tables are excluded. The raw checksum list is explicitly permitted metadata.
These commands stage the current final versions. They group the completed work for review;
they do not recreate the historical checkpoint snapshots. Test counts below describe the
actual SLURM runs. The final combined tree is verified; intermediate commits assembled now
have not been independently tested.

## 1. Real cohort, frozen split and training EDA

```bash
git add STATUS.md config/default.yaml docs/analysis-decisions.md \
  scripts/01_prepare_data.py scripts/audit_data.py scripts/eda.py scripts/slurm/run_stage.sh \
  src/brca/pipeline.py src/brca/data/load.py src/brca/data/splits.py \
  src/brca/data/prepare.py src/brca/data/preprocessing.py src/brca/data/partitions.py \
  src/brca/data/eda.py src/brca/viz/eda.py \
  tests/test_data_preparation.py tests/test_partition_access.py
git commit -m "feat(data): prepare METABRIC cohort and training EDA"
git push -u origin feat/thesis-pipeline
```

PR title: **feat(data): prepare METABRIC cohort and training EDA**

PR body: Prepare the real METABRIC cohort with DSS primary and OS sensitivity. Freeze shared
partitions, fit preprocessing on training patients only and export training EDA with source
provenance. SLURM 1185–1186 validated preparation and split reuse; independent data/leakage
review completed. Historical validation for this checkpoint: 149 passed.

## 2. Model matrix and Bayesian uncertainty

```bash
git add STATUS.md config/default.yaml pyproject.toml uv.lock src/brca/manifest.py \
  src/brca/fitting.py src/brca/models/base.py src/brca/models/cox.py \
  src/brca/models/rsf.py src/brca/models/bayesian.py \
  src/brca/evaluation/metrics.py src/brca/uncertainty.py \
  src/brca/viz/results.py src/brca/viz/uncertainty.py \
  scripts/03_fit_models.py scripts/training_uncertainty.py \
  tests/test_models_classical.py tests/test_models_bayesian.py \
  tests/test_evaluation_metrics.py tests/test_uncertainty.py tests/test_viz_results.py \
  docs/classical-models.md docs/bayesian-workflow.md docs/evaluation-implementation.md
git commit -m "feat(models): fit and diagnose the survival model matrix"
git push -u origin feat/thesis-pipeline
```

PR title: **feat(models): fit and diagnose the survival model matrix**

PR body: Implement Cox PH checks, untuned RSF and censored Bayesian AFT models for both
feature sets and endpoints. Add prior/family sensitivity, hypothetical-profile uncertainty
and verified persistence. SLURM 1204 completed 16 fits; 1210 passed 168 tests. Eight Bayesian fits
pass unchanged diagnostics. A failed numerical initialization was diagnosed, archived and
corrected with seeded finite-gradient-checked starts.

## 3. Censoring-aware conformal calibration

```bash
git add STATUS.md src/brca/conformal/survival.py src/brca/calibration.py \
  scripts/04_conformal.py tests/test_conformal.py docs/conformal-methods.md
git commit -m "feat(conformal): calibrate censored survival prediction intervals"
git push -u origin feat/thesis-pipeline
```

PR title: **feat(conformal): calibrate censored survival prediction intervals**

PR body: Add conservative restricted-time intervals and separate approximate IPCW
sensitivity, corrected ranks, explicit unsupported tails and stratified paired bootstrap.
Independent mathematical review completed. SLURM 1210 calibrated 144 settings; 1211 passed all
170 tests before final evaluation.

## 4. Final evaluation and exploratory Qin comparator

```bash
git add STATUS.md scripts/05_evaluate.py tests/test_pipeline_stages.py \
  config/qin_exploratory.yaml src/brca/conformal/qin.py tests/test_qin.py \
  docs/qin-exploratory.md config/uncertainty_reliability.yaml \
  src/brca/conformal/lower.py src/brca/conformal/qin_ridge.py src/brca/reliability.py \
  tests/test_conformal_lower.py tests/test_qin_ridge.py tests/test_reliability.py \
  docs/uncertainty-reliability.md
git commit -m "feat(evaluation): evaluate frozen and exploratory Qin survival intervals"
git push -u origin feat/thesis-pipeline
```

PR title: **feat(evaluation): evaluate frozen and exploratory Qin survival intervals**

PR body: Freeze the evaluation plan and verify caches against source, model, preprocessing
and test hashes. Report predictive metrics, paired feature comparisons, observable coverage
bounds and widths by endpoint and subgroup. SLURM 1211 passed 170 tests and completed final
evaluation. All 1,000 predictive-performance bootstrap replicates succeeded. Independent
statistical and provenance audits found no invalidating errors.
Add a separate post hoc Qin two-sided Cox bootstrap CPI for general right censoring,
reusing the verified test cache without changing registered outputs. SLURM 1332 passed
176 fast and 179 total tests; 1333 evaluated 24 exploratory interval settings. DSS tail
IPCW is unsupported and the inverse-weighted event pool is small; the comparator is not
promoted to the primary analysis.

The reliability extension passed SLURM 1337 (`make check`: 228 tests; `make test-all`:
232 tests) and completed real-data job 1338 (96 settings, all 4,000 ridge refits).
Independent audit 1341 passed 2,523 checks and final `make check` passed 228 tests.

## 5. Thesis outputs and final handoff

```bash
git add AGENTS.md Makefile README.md STATUS.md \
  docs/analysis-decisions.md docs/experimental-design.md docs/figure-style.md \
  docs/pipeline-runbook.md docs/git-checkpoints.md docs/final-review.md \
  docs/final-provenance-review.md docs/model-interface.md docs/survival-metrics.md \
  scripts/06_export_results.py scripts/audit_results.py scripts/slurm/run_pipeline.sh \
  scripts/export_reliability.py scripts/audit_reliability.py src/brca/viz/reliability.py \
  scripts/slurm/run_stage.sh src/brca/reporting.py \
  src/brca/viz/coverage.py src/brca/viz/eda.py src/brca/viz/results.py \
  tests/test_coverage_viz.py tests/test_result_audit.py tests/test_viz_results.py \
  data/raw/CHECKSUMS.sha256
git commit -m "feat(reporting): export and audit the complete thesis analysis"
git push -u origin feat/thesis-pipeline
```

Final cumulative PR title: **feat: complete METABRIC survival and uncertainty reliability analysis**

Final cumulative PR body:

Implement the real-data thesis pipeline from cohort preparation and training EDA through
Cox, Bayesian AFT and untuned RSF fitting, conformal calibration, held-out evaluation and
thesis exports. DSS is primary and OS sensitivity; all feature arms share one frozen split.
The post hoc Qin Cox bootstrap CPI is isolated from the registered primary analysis and
uses the same frozen test patients for descriptive comparison.

The analysis includes 16 fits, 144 calibration states, 168 interval-report cells, paired 1,000-draw
bootstrap comparisons and hypothetical-profile Bayesian uncertainty. All eight Bayesian
fits pass unchanged MCMC gates. Exports include 137 PDF/PNG pairs and 94 LaTeX tables.
Final SLURM job 1231 passed 171 fast tests, all 173 tests and 906 artifact integrity checks.
Independent statistical, leakage, provenance and artifact reviews are recorded in
STATUS.md and docs/final-review.md.
An additional SLURM job 1332 passed 176 fast and 179 total tests, including censored
synthetic Qin coverage, and job 1333 completed all four Qin arms. Its separate aggregate
report is documented in docs/qin-exploratory.md.

Extend uncertainty analysis with conservative one-sided lower bounds and fixed-ridge
Qin bootstraps. SLURM 1337 passed 228 fast/232 total tests, 1338 completed 96 new
interval settings and all 4,000 refits, and 1341 passed 2,523 integrity checks plus
final tests and exports. Original results are preserved. Lower-bound widths improve
in some settings, but coverage/width tradeoffs remain mixed and comparisons are post hoc.

PAM50 adds modest DSS discrimination for RSF and Weibull; Cox improvement is uncertain.
Conservative intervals have favorable marginal realised-coverage bounds but very wide
prediction intervals, and subgroup behavior is less reassuring. PH departures, net-DSS
assumptions, unsupported tail IPCW and exploratory comparisons are reported explicitly.

Generated patient/model artifacts remain excluded from Git. Dirty runs retain verified
source/config archives; clean-commit reproduction and LaTeX compilation are not claimed.

## Uncertainty reliability checkpoint after the earlier pipeline commits

If using the fresh sequence above, the implementation is already included in step 4
and its exports in step 5. If those earlier commits already exist, use this checkpoint:

```bash
git add STATUS.md config/uncertainty_reliability.yaml \
  src/brca/conformal/lower.py src/brca/conformal/qin_ridge.py src/brca/reliability.py \
  src/brca/viz/reliability.py scripts/05_evaluate.py scripts/export_reliability.py \
  scripts/audit_reliability.py tests/test_conformal_lower.py tests/test_qin_ridge.py \
  tests/test_reliability.py docs/uncertainty-reliability.md docs/conformal-methods.md \
  docs/experimental-design.md docs/git-checkpoints.md
git commit -m "feat(conformal): add survival uncertainty reliability comparisons"
git push -u origin feat/thesis-pipeline
```

PR title: **feat(conformal): add survival uncertainty reliability comparisons**

PR body: Add conservative one-sided survival bounds for all three model families and
a fixed-ridge Qin sensitivity requiring every bootstrap refit to succeed. Diagnose
censoring support and unpenalized bootstrap failures using training data only. Freeze
a separate exploratory protocol before test-cache reuse, preserve the original results,
and report paired coverage-width comparisons with subgroup summaries. Include full
comparison tables, thesis figures and an independent provenance audit. The extension
addresses numerical reliability and a different one-sided prediction question; it does
not establish clinical validity or turn the reused test set into independent validation.

Validation: SLURM 1337 passed 228 fast and 232 total tests; real job 1338 completed
all 96 settings and 4,000/4,000 ridge refits in 5m57s. Final job 1341 passed 228 fast
tests and 2,523 independent source/artifact checks, then exported the report and four
PDF/PNG pairs. All four previews were visually inspected. Original plans and results
retain their historical manifest hashes.

## Standalone Qin checkpoint if steps 1–5 were already committed

If the earlier five commits already exist, stage only this later comparator and its
documentation. If using the fresh sequence above, Qin is already included in step 4
and this standalone command should be skipped.

```bash
git add STATUS.md docs/conformal-methods.md docs/git-checkpoints.md \
  docs/qin-exploratory.md config/qin_exploratory.yaml \
  src/brca/conformal/qin.py scripts/05_evaluate.py tests/test_qin.py
git commit -m "feat(conformal): add exploratory Qin bootstrap survival intervals"
git push -u origin feat/thesis-pipeline
```

PR title: **feat(conformal): add exploratory Qin bootstrap survival intervals**

PR body: Implement Qin et al.'s two-sided bootstrap CPI with a Cox working model,
inverse-censoring resampling of observed training events and a separate, source-bound
evaluation using the verified original test cache. Preserve the registered primary
plan and metrics. SLURM 1332 passed 176 fast and 179 total tests; job 1333 completed
four arms and 24 interval settings. Report coverage bounds with interval widths,
and qualify the weak DSS censoring support, Cox refit failures and post hoc test reuse.
