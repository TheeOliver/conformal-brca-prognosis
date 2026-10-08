# Final thesis-analysis review — 2026-09-29

The implemented analysis uses the real METABRIC clinical extract with 1,979 eligible patients,
DSS primary and OS sensitivity. The one frozen split contains 989/495/495 patients. Models,
preprocessing, metric grids, target horizons and calibration were fixed before test access.
The final report is `outputs/reports/thesis_results.md`; machine-readable RQ evidence is
`outputs/metrics/research_answers.json`. Generated artifacts remain outside Git.

## Verification record

| Check | Evidence |
| --- | --- |
| Source schema, cohort, preprocessing, split/EDA isolation | Independent review; SLURM 1185/1186 |
| Real model matrix | SLURM 1204: 12 primary fits and 4 training-only sensitivities completed |
| Bayesian acceptance | All 8 fits: zero divergences; worst R-hat 1.0059; minimum bulk ESS 1714, tail ESS 1956 |
| Calibration | SLURM 1210: 144 states; 168 fast tests passed; training-profile outputs generated |
| Full final evaluation | SLURM 1211: 170 tests passed, 495 test patients, 1,000 paired bootstrap replicates |
| Mathematical/caller review | Independent score domination, corrected ranks, censoring ties, primary-stratum and no-leakage checks |
| Aggregate statistical review | SLURM 1217/1219: no conclusion-invalidating findings |
| Export and report generation | SLURM 1224: 168 tests passed; 378 export artifacts recorded |
| Artifact inspection | SLURM 1226: 137 PDF/PNG pairs, all PDFs 6.3 in wide; 94 booktabs tables structurally checked |
| Final independent integrity audit | SLURM 1230: 906 checks, zero failures; see `final-provenance-review.md` |
| Integrated final verification | SLURM 1231: `make format`, `make check` (171 passed), `make test-all` (173 passed), provenance audit (906 checks, zero failures) |

Independent reviewers worked in separate worktrees. Findings were corrected before this
handoff: bootstrap failure reporting and sampling strata; source/cache bindings; explicit
PH caveats; exact metric horizon; unsupported-tail reporting; precise paired IBS intervals;
figure dimensions, units, axis labels and grayscale distinctions. No models or calibration
were changed after test access. Subsequent edits affected reporting and visualization only.

One initial Bayesian OS molecular chain was numerically stuck despite a finite initial log
probability. Its gradient was non-finite. The failed posterior was preserved and diagnosed;
bounded seeded starts with finite-gradient checks resolved it. All 16 fits were then rerun
under common settings, without changing the priors or loosening diagnostic thresholds.

## Answers and limits

1. **Clinical prediction:** DSS Uno C is approximately 0.686–0.690 across the three clinical
   models; the report contains bootstrap intervals and accompanying IBS/Brier/AUC.
2. **PAM50:** paired DSS Uno C gains favor RSF and Weibull; Cox's interval includes zero.
   RSF IBS improves modestly; Weibull IBS improvement remains uncertain. These are exploratory
   internal-validation comparisons without multiplicity adjustment.
3. **Model families:** Bayesian molecular and RSF molecular have similar Uno C point estimates
   near 0.705; the report does not claim family superiority. All four Cox PH diagnostics flag
   departures. Weibull AFT also imposes PH; the lognormal training sensitivity is retained.
4. **Bayesian uncertainty:** posterior curve bands and individual predictive time intervals
   answer different questions. Family/prior sensitivity and PPCs matter despite passing MCMC.
   Hypothetical profiles may have weak joint covariate support; unrestricted tails are model
   extrapolations, not literal clinical forecasts.
5. **Conformal coverage:** primary conservative marginal realised-coverage lower bounds exceed
   nominal in this sample at all three levels. At 300 months, median widths span 251–287 months.
   DSS tail IPCW is unavailable because estimated censoring survival is 0.01478, below 0.05.
   Marginal behavior does not guarantee subgroup behavior: several Basal subgroup upper bounds
   lie below nominal, with wide uncertainty intervals.
6. **Discrimination versus uncertainty:** predictive ranking alone does not establish useful,
   calibrated individual uncertainty. Wide intervals and subgroup variability remain material
   limitations. The report includes a joint C-index/coverage/width table at every nominal level.

IBS integrates over 28.973–120 months in this run. Bootstrap intervals preserve primary event
strata and condition on fixed fits, calibration and censoring estimation; they omit training
and split-selection variability. Net DSS requires explicit competing-death/censoring
interpretation. External validation remains future work, not a result of this project.

## Reproducibility and remaining integration checks

Runs explicitly record dirty Git state and retain content-addressed source/config archives,
config/artifact hashes, resolved package versions, SLURM jobs and elapsed time. Review checked
those archived bytes; it did not claim reproduction from a clean Git commit. User checkpoint
commands and PR text are in `git-checkpoints.md`. The raw checksum list is safe metadata;
patient tables, caches and posterior files must never be staged.

PDF/PNG dimensions, readability and table structure were checked. No LaTeX engine is installed
on the cluster, so compilation within the external thesis document remains unverified.
Tables require `booktabs`, `siunitx`, and `longtable` for coverage appendices. The existing
Codex-hook runtime integration remains unverified; in-code compute guards and SLURM execution
were exercised and do not depend on those optional hooks.
