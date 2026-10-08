# Survival-uncertainty reliability experiment

This post hoc experiment addresses two weaknesses exposed by the original results:
long, weakly supported survival intervals and failed unpenalized Cox bootstrap refits.
Its objective is more dependable evidence about uncertainty, not selecting a new primary
method by whichever test result looks best. The primary analysis and first Qin comparator
remain frozen. All candidates use the original 989/495/495 DSS-event-stratified split,
both endpoints, both feature sets, miscoverage 0.20/0.10/0.05, and horizons 120/300 months.

## One-sided conservative survival bounds

For a fixed fitted model, let `m(x)` be its alpha quantile, clipped at the horizon `tau`.
The signed calibration score is `m(X) - min(Y, tau)`, with observed follow-up
`Y = min(T, C)`. It is at least the unobserved score `m(X) - min(T, tau)` for every
patient, including censored patients. Taking the exact order statistic at
`ceil((n + 1) * (1 - alpha))` therefore gives a conservative lower prediction bound.
An unavailable rank returns full support. Signed negative thresholds are retained.

Because the frozen split stratifies on the primary DSS event indicator, calibration
computes the threshold within each of those strata and takes the larger threshold.
The model does not need the future patient's stratum. This gives a valid score-envelope
argument under exchangeability within each split stratum for a rule fixed before
calibration. It does not require independent censoring for the score inequality.
It does not provide conditional coverage for an individual patient, simultaneous
coverage across all analyses, or protection after test-driven method selection.

The resulting set `[L(x), tau]` targets `min(T, tau)`. Equivalently, its promise is
`T >= L(x)`: survival of at least this long. `tau` is an administrative target boundary,
not a predicted upper survival limit. Its width `tau - L(x)` describes remaining
restricted support. The new method changes the prediction question; a narrower set
does not by itself prove superiority to a two-sided Qin interval.

All three existing model families supply the base quantile, without refitting.
Only the 495 calibration patients determine the lower-bound correction. The method
is a conservative score-envelope adaptation, not the published Candes/Gui algorithm.
For DSS, the latent net-survival target still requires careful interpretation when
other-cause deaths preclude an observed breast-cancer death.

## Fixed-ridge Qin sensitivity

The [original Qin comparator](qin-exploratory.md) produced failed bootstrap Cox fits,
especially with molecular features. The new sensitivity fixes the Cox ridge penalty
at **1.0**, using the existing training-only preprocessing. This choice was fixed before
this experiment accessed the already-opened test cache; it was not tuned using test
coverage. Original and bootstrap estimators use identical penalties and solver settings.
All **1,000 requested refits must succeed**. Any failure rejects the complete candidate;
there is no distribution conditioned on successful fits.

The resampling and censoring weights otherwise follow the original implementation of
[Qin et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12104816/).
Identical named seeds preserve paired training resamples and sampled events. A separate
100-resample unpenalized training pilot records constant columns, rank deficiencies and
exception types. It diagnoses numerical problems and is not a new calibration sample.

Ridge is an exploratory working-model adaptation. It supplies no new coverage theorem,
does not repair PH departures, and does not increase the effective size of the
inverse-censoring-weighted event pool. The demanding marginal censoring assumption
remains. The training support audit reports those limitations separately at each horizon.

## Evaluation and reproducibility

Before preparation, a separate immutable plan binds original inputs, source, lockfile,
settings, and reference results. All fitting and calibration finish before stage 05
reads the hash-verified original test cache. Prepared models and thresholds have their
own binding. Reconstructed historical interval statistics must match the saved originals
before paired comparisons are accepted.

The experiment reports 96 interval settings: 72 lower-bound settings and 24 ridge-Qin
settings, with coverage bounds, width, supported IPCW estimates, subgroup summaries,
and paired 1,000-draw test bootstrap differences. Bootstrap intervals condition on fitted
models, calibration and training censoring estimates; they do not include refitting
uncertainty. Subgroups smaller than 10 are suppressed. Multiple comparisons are
exploratory. The original test set has already informed development, so fresh external
validation is needed before selecting a new primary procedure.

Ridge predictive metrics use the original training-selected evaluation grids and shared
test bootstrap protocol. Feature-set differences are paired; ridge-versus-original
predictive point comparisons do not have a dedicated paired significance test.

```bash
sbatch scripts/slurm/run_stage.sh 05 config/default.yaml --allow-dirty --reliability-exploratory
```

After successful evaluation, export the saved aggregates and verify their bindings:

```bash
sbatch --partition=openlab-queue --cpus-per-task=2 --mem=8G --time=00:10:00 \
  --output=logs/reliability-export-%j.out \
  --wrap='uv run python scripts/export_reliability.py --allow-dirty'
sbatch --partition=openlab-queue --cpus-per-task=1 --mem=2G --time=00:10:00 \
  --output=logs/reliability-audit-%j.out \
  --wrap='uv run python scripts/audit_reliability.py --project-root "$PWD" --job 1338 --output outputs/uncertainty_reliability/independent_audit.json'
```

The audit command above names this experiment's job; use the corresponding evaluation
job ID when auditing a separately reproduced run. The exporter generates
`outputs/uncertainty_reliability/comparison_report.md` and four PDF/PNG pairs in
`outputs/figures/`. It uses saved aggregate statistics only. The independent audit
checks archived sources, original manifest hashes, complete experiment cells, bootstrap
accounting, subgroup suppression and paired point differences; it hashes patient
artifacts without deserializing their contents.

Outputs are isolated under `outputs/uncertainty_reliability/`. No patient-level artifact
is versioned. Changed statistical source deliberately prevents the original opened
analysis from being rerun in place; its archived source and frozen results remain the
authoritative original analysis. Clean-commit numerical reproduction remains separate.

## Results and achieved improvement (2026-10-08)

SLURM **1338** completed all 96 interval settings in **5 minutes 57 seconds**.
The independent audit in job **1341** passed **2,523/2,523** source, artifact and
aggregate-consistency checks, including preservation of the original primary and Qin
plans/results. This is an integrity audit, not a second numerical reproduction.

### Qin fits are now numerically dependable

| Endpoint | Features | Original successful refits | Ridge successful refits |
| --- | --- | ---: | ---: |
| DSS | Clinical | 1,000/1,000 | 1,000/1,000 |
| DSS | Clinical + PAM50 | 900/1,000 | 1,000/1,000 |
| OS | Clinical | 998/1,000 | 1,000/1,000 |
| OS | Clinical + PAM50 | 931/1,000 | 1,000/1,000 |

All **4,000** planned ridge refits succeeded, recovering the 171 draws lost by the
original fits. The unpenalized 100-draw diagnostic pilot had 9 DSS-molecular and 11
OS-molecular failures; respectively 4 and 8 coincided with a constant/rank-deficient
design. The remaining failures raised numerical warnings. The ridge runs retained
44 DSS-molecular and 47 OS-molecular resamples with absent unknown-PAM50 categories
and completed them successfully. This supports numerical stabilization, not a claim
that absent categories explain every original failure.

At 90% nominal coverage, ridge Qin median widths differ from original Qin by less
than one month across all eight endpoint/feature/horizon settings. For DSS + PAM50
at 120 months, width **increases** from 91.1 to 92.0 months, while observable coverage
bounds change from 91.7–94.1% to 92.1–94.5%. The gain is a complete bootstrap
distribution rather than a dramatic increase in predictive precision.

### Lower bounds offer a useful additional question, with mixed tradeoffs

The table below uses DSS, **90% nominal coverage**, the **120-month horizon** and
clinical + PAM50 predictors. Each lower-bound calibration uses 495 patients and is
evaluated on the same 495 test patients. Widths are months; confidence intervals are
95% paired test-bootstrap intervals conditional on the fixed fits and calibrations.

| Base model | Primary median width | Lower-bound median width | Paired width change [CI] | Lower-bound observable coverage range | Lower-bound IPCW coverage [CI] |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cox | 94.0 | 88.1 | -5.9 [-6.3, -4.6] | 88.9–94.7% | 94.4% [92.4, 96.3] |
| Weibull | 95.1 | 86.4 | -8.8 [-9.9, -8.1] | 89.5–95.2% | 94.8% [92.8, 96.8] |
| RSF | 91.4 | 87.6 | -3.8 [-5.8, 0.3] | 87.5–95.2% | 94.8% [92.8, 96.7] |

These are narrower restricted prediction sets, not demonstrated gains in two-sided
prediction. The observable ranges cross 90%, so censoring alone does not establish
that their realised latent test coverage reaches nominal. IPCW estimates above nominal
depend on the censoring assumptions. The RSF width-change interval includes zero.
At 90% nominal coverage, DSS lower-bound median widths across both feature arms are
85.5–90.1 months at the 120-month horizon, and **265.5–270.1 months** at 300 months:
still broad. OS is mixed: the clinical Cox lower-bound width increases by 4.8 months
at 120 months; clinical Weibull increases by 2.7 months, whereas RSF decreases by
5.3–5.6 months across feature arms. Other nominal levels are all retained in the report.

Marginal behavior does not imply reliable PAM50 subgroup coverage. For example, the
DSS Cox + PAM50 lower bound at 90% nominal/120 months has an observable coverage range
of **84.7–86.4%** among 59 Basal patients, despite median width **101.8 months**. The
observable upper bound's bootstrap interval is 79.7–94.9%, so this identifies low
realised subgroup coverage without establishing population subgroup undercoverage.
All prespecified subgroup summaries remain available; tiny groups are suppressed.

### Tail support remains the central limitation

DSS training censoring survival is **0.672 at 120 months**, but **0.0148 at 300 months**.
Only 7 training patients have observed follow-up through 300 months, compared with
471 through 120 months. The registered 0.05 support gate therefore still suppresses
300-month DSS IPCW coverage. The full Qin DSS event pool retains effective size **12.8**
despite 323 observed training events; ridge does not alter those weights.

The project now has a stable Qin sensitivity, a separately justified lower-bound method,
and evidence showing where precision and identifiability remain limited. There is no
universal winning method, no replacement of the primary analysis, and no new claim of
individual clinical reliability. All new comparisons remain post hoc.

### Verification and handoff

- SLURM **1337**: `make format`, `make check` (**228 passed**) and `make test-all`
  (**232 passed**), including synthetic censoring/coverage tests and full synthetic
  stage integration with provenance-tampering checks.
- SLURM **1341**: final `make check` (**228 passed**), independent audit
  (**2,523 passed**) and saved-aggregate exports (full Markdown report, four PDF/PNG
  pairs, export metadata and manifest). All four PNG previews were visually inspected.
- Independent mathematical, leakage and numerical-method reviews found no remaining
  confirmed implementation defect. Synthetic coverage checks are regression evidence;
  the ridge test covers only three generated populations and is not broad validation.
- Intermediate jobs 1334–1336 failed on test-fixture lint/import/manifest-lifecycle
  issues; 1339–1340 stopped on audit/export lint. These were corrected before the
  successful checks above. The completed real experiment was not rerun or modified.
- No patient data or generated results were staged; no commit, push or PR was created.
  Exact review/commit commands and PR text are in [git-checkpoints.md](git-checkpoints.md).

Full aggregate report: [comparison_report.md](../outputs/uncertainty_reliability/comparison_report.md).
Independent audit: [independent_audit.json](../outputs/uncertainty_reliability/independent_audit.json).
