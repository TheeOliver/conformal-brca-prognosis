# Qin et al. bootstrap CPI: exploratory comparison

This is a **post hoc comparator** to the frozen primary METABRIC analysis, not a replacement
for its prespecified split-conformal score envelope. It uses the same training and test
patients, endpoints, feature sets, alpha values, and restricted horizons. The saved original
Cox models and training-only preprocessors are reused; no model is chosen using the test set.
The original calibration patients are not used because Qin et al.'s resampling algorithm
operates on the training data. The original test cache is reused after its SHA-256 check.
Consequently, differences seen on test are descriptive and cannot validate a post hoc
method selection.

The implemented two-sided **Cox working-model** variant follows Qin, Piao, Ning and Shen,
“Conformal predictive intervals in survival analysis: a resampling approach,” *Biometrics*
(2025), DOI [10.1093/biomtc/ujaf063](https://doi.org/10.1093/biomtc/ujaf063).
Their [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12104816/) and
[author implementation](https://github.com/JPiao7u089/CPI-bootstrap) describe the algorithm.
The Python implementation here was written independently. It does **not** implement the
paper's AFT variants.

For each endpoint and feature set, estimate the marginal reverse-Kaplan–Meier censoring
survival \(G(t-)\) on training patients. Give each observed training event a sampling mass
proportional to \(1/G(Y_i-)\). On each bootstrap iteration, resample all training patients,
refit the Cox model, sample one event patient using those masses, and record the bootstrap
survival pivot \(S_b(Y_i\mid X_i)\). The lower and upper pivot quantiles are the ordinary
sample \(\alpha/2\) and \(1-\alpha/2\) quantiles used by the paper, **not** the finite-rank
split-conformal correction. Invert the original fitted Cox survival curve at those two
thresholds to obtain a time interval. Clip both endpoints at each registered horizon to
target \(\min(T,\tau)\). Cox baseline survival after its last event knot is carried forward,
as in the author's step-function implementation; the number of bootstrap pivots needing
this tail convention is reported. The upper time edge is the first knot where survival is
**strictly below** the lower pivot threshold, so an equality plateau stays in the set;
reporting the closed endpoint is a conservative closure at that knot. For an event and
censoring tie, the inverse-censoring calculation uses \(G(Y-)\), treating the event as
observed before censoring at that time. This explicit tie convention can differ from the
author's R `survfit` evaluation on discretely tied clinical months.

This method needs only general right-censored observations \((Y,\delta,X)\); it does not
need every patient's unobserved censoring time. Its validity is **approximate**, subject
to censoring assumptions and adequate support. Our marginal \(G\) assumes censoring is
independent of both event time and covariates, which is demanding for disease-specific
survival when deaths from other causes are treated as censoring. The saved preprocessor is
held fixed across bootstrap refits, an explicit adaptation of the author's raw-covariate
procedure. Cox proportional-hazards diagnostics already flag departures in the original
analysis. No finite-sample or individual-patient coverage guarantee is claimed.

All outcomes on the test set remain partially censored. The report therefore gives
observable lower/upper bounds on realised coverage, an IPCW point estimate only where the
registered censoring-support gate permits it, and interval-width quantiles. The same
test-patient bootstrap protocol is used for uncertainty as in the primary report. At long
disease-specific follow-up, IPCW coverage may be unavailable. `n_cal` in this exploratory
report is the number of successful **bootstrap pivots**, not held-out calibration patients.

Run only through SLURM:

```bash
sbatch scripts/slurm/run_stage.sh 05 config/default.yaml --allow-dirty --qin-exploratory
```

The isolated output is `outputs/qin_exploratory/coverage.json` and its plan is
`outputs/qin_exploratory/plan.json`. Neither overwrites `outputs/evaluation/plan.json`
or any registered primary metric. The changed stage-05 source means rerunning the original
stage-05 script against its opened plan will fail its source-hash check; its archived
original source and frozen outputs remain the authoritative original run.

## METABRIC run and interpretation (2026-10-08)

SLURM job **1333** completed all four Cox arms with 1,000 requested bootstrap refits
each. The independent read-only implementation review found no remaining confirmed
algorithm or leakage defect after corrections. SLURM job **1332** passed `make check`
(176 fast tests) and `make test-all` (179 total, including a censored synthetic
end-to-end coverage test). The run manifest is `outputs/manifests/05_qin_exploratory_*.json`.
The earlier job 1329 was canceled during review; its plan is preserved as
`outputs/qin_exploratory/plan_cancelled_1329.json` and contributes no reported result.

The table gives **90% nominal** intervals on the same 495 held-out patients in each
endpoint. “Observable range” bounds realised latent coverage; it is not a confidence
interval. Width is the median in months. All 80%, 90%, and 95% marginal and subgroup
results, including fixed-fit test bootstrap intervals, are in the aggregate JSON.

| Endpoint | Features | Horizon (months) | Qin observable range | Qin IPCW | Qin width | Primary conservative Cox width |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| DSS | Clinical | 120 | 92.3–94.7% | 94.4% | 92.7 | 94.5 |
| DSS | Clinical + PAM50 | 120 | 91.7–94.1% | 93.6% | 91.1 | 94.0 |
| DSS | Clinical | 300 | 90.5–94.7% | unavailable | 272.1 | 274.4 |
| DSS | Clinical + PAM50 | 300 | 89.3–94.1% | unavailable | 270.3 | 273.1 |
| OS | Clinical | 120 | 96.0–97.0% | 96.7% | 101.0 | 97.4 |
| OS | Clinical + PAM50 | 120 | 91.3–92.7% | 92.2% | 93.7 | 99.8 |
| OS | Clinical | 300 | 89.5–97.0% | 96.5% | 277.8 | 274.8 |
| OS | Clinical + PAM50 | 300 | 85.3–92.3% | 90.9% | 270.3 | 277.2 |

The Qin intervals are only modestly shorter than the conservative primary intervals
at the DSS 300-month target. The 300-month DSS IPCW estimate is unavailable under the
registered 0.05 censoring-support gate. More seriously, the smallest *training-event*
censoring survival is about **0.0049**, so the 323 observed DSS events have an
inverse-weighted effective event-pool size of only **12.8**. The marginal censoring
assumption is also questionable for DSS because other-cause deaths count as censoring.
Thus these are descriptive computed intervals, not evidence that Qin's asymptotic
coverage guarantee applies well to this cohort or tail.

Bootstrap Cox refits succeeded in 1,000/1,000 DSS clinical, **900/1,000 DSS
clinical+PAM50**, 998/1,000 OS clinical, and 931/1,000 OS clinical+PAM50. The
DSS molecular arm is exactly at the prespecified 90% success gate; its successful-fit
conditioned pivot distribution deserves particular caution. Original PH diagnostic
departures and near-full-horizon widths further limit clinical usefulness. These
comparisons reuse the test set after the original primary analysis and cannot justify
choosing a new primary method by which coverage looks better.

The registered original evaluation plan and coverage JSON retained their original
manifest SHA-256 hashes. The Qin run used an explicitly dirty working tree, with its
exact code/config captured in a source archive; clean-commit numerical reproduction
remains outstanding.
