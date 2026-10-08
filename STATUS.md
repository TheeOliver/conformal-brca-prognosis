# STATUS

Shared progress ledger for **Claude Code and Codex**. Read it at the start of every session;
update it at the end of every session that changed files. Procedure:
`.claude/skills/status-handoff/SKILL.md`. Newest handoff entries go at the **top** of the log.

## Snapshot

- **Phase:** 6 — primary analysis complete; post hoc uncertainty reliability extension audited
- **Last updated:** 2026-10-08 by Codex
- **Next up:** user checkpoint commits/PR (`docs/git-checkpoints.md`), supervisor review of
  `docs/uncertainty-reliability.md`, thesis integration and clean-commit numerical reproduction.

## Active claims

Claim before starting non-trivial work; remove your row when you hand off.

| Task | Agent | Branch | Since |
| --- | --- | --- | --- |

## Milestones

**Phase 0 — Setup**
- [x] Agent configuration: `AGENTS.md` (canonical), `CLAUDE.md` (imports it), rules, skills,
      review checklists, hooks for Claude and Codex
- [x] Agent-independent guards: `.githooks/pre-commit`, `src/brca/compute_guard.py`
- [x] Environment: uv, Python 3.12, `uv.lock`
- [x] Synthetic fixture + scaffold tests (`make check` → 19 passed)
- [x] First commit — *user* (`d22268b`, now on `main` / `origin/main`)

**Phase 1 — Data**
- [x] METABRIC clinical data acquired: cBioPortal public study, 1,981 usable patients
      (`make download`; `docs/metabric-data-dictionary.md`)
- [x] Endpoint decided: DSS primary; OS sensitivity (user, 2026-09-29)
- [x] Reproducible aggregate data audit before cohort decisions
- [x] Stage 01 `prepare_data` + tests
- [x] Stage 02 `make_splits` + tests (synthetic frozen `splits.json`, disjointness asserted)
- [x] Freeze the real-data split after cohort decisions
- [x] Training-set EDA: predictor distributions, missingness, censoring and survival patterns

**Phase 2 — Models** (each for *both* feature sets)
- [x] Cox PH, with proportional-hazards check
- [x] Bayesian Weibull AFT: prior predictive, diagnostics gate, PPC, prior sensitivity
- [x] Random Survival Forest (untuned, by design)

**Phase 3 — Conformal**
- [x] Censoring adaptation chosen and justified (see `docs/conformal-methods.md`)
- [x] Split CP with finite-sample quantile + coverage unit tests
- [x] `conformal-validity-auditor` review — independent mathematical and caller reviews

**Phase 4 — Evaluation**
- [x] Uno C, Brier, IBS, td-AUC with bootstrap CIs on a shared interior time grid
- [x] Coverage + width at 80/90/95%, marginal and by PAM50 / ER subgroup

**Phase 5 — Answers** (generated artifacts remain ignored)
- [x] RQ1: clinical prediction — [performance](outputs/metrics/performance.json)
- [x] RQ2: added PAM50 information — [paired comparisons](outputs/metrics/rq2_feature_set_comparison.json)
- [x] RQ3: model families — [performance](outputs/metrics/performance.json)
- [x] RQ4: Bayesian uncertainty — [answers and evidence](outputs/metrics/research_answers.json)
- [x] RQ5: coverage with width — [coverage report](outputs/metrics/conformal_coverage.json)
- [x] RQ6: discrimination versus uncertainty — [joint evidence](outputs/metrics/research_answers.json)

**Phase 6 — Thesis artefacts**
- [x] Figures and tables exported per `docs/figure-style.md`: 137 PDF/PNG pairs, 94 LaTeX tables
- [x] `leakage-auditor`, `stats-reviewer`: independent reviews completed, findings corrected
- [x] `reproducibility-checker`: 906 artifact/source integrity checks passed
- [ ] Strict clean-commit numerical reproduction sign-off: archived dirty sources verified,
      but no second real-data numerical run or clean-commit reproduction is claimed

**Post hoc published-method comparator**
- [x] Qin et al. two-sided bootstrap CPI with Cox working model, train-only refits and
      inverse-censoring event resampling; synthetic censored coverage test and independent review
- [x] Separate four-arm DSS/OS × feature-set SLURM evaluation, 80/90/95% intervals at
      120/300 months, with coverage bounds, IPCW where supported, width and subgroup reports

**Post hoc uncertainty reliability improvement**
- [x] Conservative signed-score lower bounds across all three model families, with
      finite-rank, stratification and censoring-coverage synthetic tests
- [x] Fixed-ridge Qin sensitivity: matched original/bootstrap estimators and strict
      1,000/1,000 refit requirement; training support and failure diagnostics
- [x] SLURM 1338: 96 interval settings and paired method/feature comparisons on frozen
      patients, with original reference statistics verified before comparison
- [x] SLURM 1341: 2,523 independent integrity checks; full aggregate report and four
      PDF/PNG comparison pairs; final `make check` 228 passed (232 total in job 1337)

Results: [thesis report](outputs/reports/thesis_results.md). Verification and qualifications:
[final review](docs/final-review.md), [independent provenance review](docs/final-provenance-review.md).

## Blockers & open questions

**2026-09-29 resolution:** the four stage-01 questions below are superseded by the user's
DSS/OS choice and delegated decisions in `docs/analysis-decisions.md`: retain claudin-low,
NC → unknown, train-only imputation, full follow-up. They no longer block the pipeline.

- ~~**BLOCKER — METABRIC access.**~~ Resolved 2026-09-29: cBioPortal `brca_metabric` is a
  public study with all needed clinical fields + PAM50 + OS/RFS. *(Claude)* Simona used a
  different distribution of the same study (Kaggle, 1,904 patients); her cohort definition
  reproduces on ours (same 2 exclusions, 32.6% vs 32.7% DSS event rate). See
  `docs/metabric-data-dictionary.md` § Cross-check.
- ~~**Endpoint — blocks stage 01.**~~ Resolved: DSS primary, OS sensitivity. Historical question:
  497 of 1,144 deaths (43%) are from other causes, so overall
  vs disease-specific survival will give materially different models; relapse-free is also
  available. Decide with the supervisor. **Precedent: Simona used disease-specific** (OS only
  as sensitivity). *Owner: user.*
- ~~**PAM50 levels — blocks stage 01.**~~ Resolved: retain claudin-low, NC → unknown.
  Historical question: data has claudin-low (218) and NC (< 10) beyond the 5
  schema levels. Proposal: claudin-low as its own level (**as Simona did**), NC (6) →
  `unknown`. *Owner: user.*
- ~~**Missing predictors — blocks stage 01.**~~ Resolved: training-only numeric median plus
  missingness indicators; categorical unknown. Historical question: 166 of 1,981 lack a predictor (grade 87, nodes 76,
  size 25). Complete-case or train-only imputation? Simona imputed (MICE for size, mode for
  grade), apparently before her split — ours must be train-only. *Owner: user.*
- ~~**Admin censoring.**~~ Resolved: retain full follow-up for fitting; restrict the primary
  conformal target to 300 months, with prespecified 120-month sensitivity. Historical question:
  config says 300 months; follow-up runs to 355. Simona did not truncate
  and capped IBS at 120 months. Confirm. *Owner: user.*
- ~~**AFT likelihood (Phase 2).**~~ Resolved: Weibull primary; lognormal and alternative-prior
  training sensitivities completed, with PPCs and explicit adequacy limitations. Historical question:
  Simona found Weibull the worst-fitting AFT family by AIC
  (log-normal best). The brief asks for Bayesian Weibull/AFT — compare against log-normal by
  posterior predictive checks before committing. *Owner: user + whoever builds the model.*
- ~~**Conformal censoring adaptation.**~~ Resolved: conservative restricted-time score envelope,
  with separate approximate IPCW sensitivity; mathematical and caller reviews completed.
  Historical question: Candès-style lower predictive bound, IPCW-weighted scores,
  or Gui adaptive cut-offs — see `docs/conformal-methods.md`. Decide before Phase 3.
- ~~**Ask Simona:** preprocessing, missing-PAM50 handling, endpoint.~~ Answered by her thesis
  (`references/simona.md`) except NC handling and whether her imputation was fit before the
  split — worth asking her those two.

## Decisions

| Date | By | Decision | Why |
| --- | --- | --- | --- |
| 2026-09-20 | Claude | Python 3.12 via uv; `uv.lock` committed | reproducible env on the cluster |
| 2026-09-20 | Claude | Split 50/25/25, stratified by event, seed 8927 | calibration set sized for α = 0.05 |
| 2026-09-20 | Claude | `pm.Censored` is the only censoring idiom | `pm.Potential` fails silently when wrong |
| 2026-09-20 | Claude | Uno's C (`concordance_index_ipcw`), not Harrell's | Harrell's is biased under heavy censoring |
| 2026-09-20 | Claude | RSF left untuned | brief rules out hyperparameter optimisation |
| 2026-09-20 | Claude | Pipeline stages 01–05 are stubs until implemented | analysis code is the thesis work |
| 2026-09-20 | user | Commit/push reserved for the user; no AI attribution in commits or PRs | user preference |
| 2026-09-29 | Claude | `AGENTS.md` canonical; `CLAUDE.md` imports it | one source of conventions, no drift |
| 2026-09-29 | Claude | Critical guards live in git hooks + code, not agent hooks | must hold for both agents and humans |
| 2026-09-29 | Claude | Synthetic fixture generated on first test run, never committed | keeps "no `.csv` in git" absolute |
| 2026-09-29 | Claude | Data source: cBioPortal API, clinical tables only; raw kept as gzipped JSON exactly as returned | public; no expression data needed; raw must stay untouched |
| 2026-09-29 | Claude | `data/raw` snapshot is immutable — `make download` refuses to overwrite | a silently replaced snapshot makes results untraceable |
| 2026-09-29 | Claude | The 528 patients without outcome/PAM50 (2016 sequencing extension) are out of scope | cannot enter a survival model |
| 2026-09-29 | Codex | Stage 02 implemented; stages 01 and 03–05 remain stubs (supersedes earlier blanket stub status) | stage 02 was independently unblocked |
| 2026-09-29 | Codex | Split JSON schema v1 has protocol/source metadata plus `splits` containing train/calibration/test ID lists; assignments use sorted IDs and a seeded Generator | shared partitions for both feature sets, independent of row order |
| 2026-09-29 | Codex | Exact sizes use largest-remainder rounding; event counts use bounded proportional allocation, with both classes required in each set | deterministic integer allocation, including feasible rare-event cohorts |
| 2026-09-29 | Codex | Compatible stage-02 reruns reuse without resampling; changed cohort bytes or split protocol fail without overwrite | protect downstream fits/calibration; other config changes get a new manifest hash |
| 2026-09-29 | Codex | Synthetic mode is explicit; dirty/unborn Git runs require `--allow-dirty` and record the provenance limitation | development checks are not thesis results |
| 2026-09-29 | user | DSS primary, OS sensitivity; remaining cohort choices delegated | explicit endpoint decision |
| 2026-09-29 | Codex | Retain claudin-low; NC/missing subtype → unknown; training-only median numeric imputation with indicators and categorical unknown; retain full follow-up | retain eligible patients without fitting preprocessing outside training; `docs/analysis-decisions.md` |
| 2026-09-29 | Codex | All stages 01–06 implemented and run on real data; supersedes prior stub decisions | 1,979 eligible patients, frozen 989/495/495 split |
| 2026-09-29 | Codex | 12 main fits across endpoints/feature arms plus 4 DSS Bayesian sensitivities; no tuning on test | registered common comparison and training-only adequacy checks |
| 2026-09-29 | Codex | Seeded bounded prior-centered MCMC starts with finite-gradient checks; gates and priors unchanged | resolved a diagnosed numerical initialization failure before test access |
| 2026-09-29 | Codex | Conservative CP targets min(T, 300 months), with 120-month sensitivity; finite ranks within primary DSS strata and worst-stratum threshold | censoring score envelope preserves conservative comparison under stated exchangeability assumptions |
| 2026-09-29 | Codex | IPCW is a separate approximation and unavailable below training censoring support 0.05; always report observable bounds and width | avoid asserted latent coverage in unsupported tails |
| 2026-09-29 | Codex | Freeze source/config/model/calibration plan before test access; paired 1,000-draw bootstrap preserves primary DSS strata | no post-test changes to fitting/calibration; fixed-fit uncertainty is explicitly qualified |
| 2026-09-29 | Codex | Actual IBS integration grid is 28.973–120 months, selected using training data | report the evaluated range exactly |
| 2026-09-29 | Codex | Dirty runs retain immutable source/config/lock archives; integrity verification is distinct from clean-commit numerical reproduction | user reserves commits; final audit cannot claim an unperformed reproduction run |
| 2026-10-08 | Codex | Qin's two-sided Cox bootstrap CPI is post hoc and isolated under `outputs/qin_exploratory/`; it reuses frozen train/Cox fits and verified test cache, not the calibration patients | original primary plan/results stay authoritative; test reuse cannot select a new primary method |
| 2026-10-08 | Codex | Qin uses 1,000 requested refits with a 90% success gate, marginal reverse-KM `G(t-)` event weights and fixed train-only preprocessing | paper's general right-censoring algorithm with explicit finite-sample and censoring-assumption limitations |
| 2026-10-08 | Codex | Separate exploratory reliability plan fixes ridge alpha 1.0, 1,000/1,000 required refits, 100-resample unpenalized diagnostic pilot; original Qin unchanged | remove conditioning on successful refits while preserving the original comparator and its limitations |
| 2026-10-08 | Codex | Add signed-score conservative lower bounds for all 12 frozen fits, with worst-primary-event-stratum thresholds; same alpha/horizon grid | answer T ≥ L under the restricted target without unstable censoring weights; this is a different one-sided question |
| 2026-10-08 | Codex | All new results remain post hoc; preserve original plans, fits and metrics; immutable source/input plan and preparation binding precede additional test-cache read | reuse cannot validate test-driven selection; no replacement of primary 300-month DSS analysis |

## Known issues

- **Codex hooks unverified.** `.codex/hooks.json` targets Codex 0.153.4's hook format, but
  could not be exercised from this session. 30-second probe: open Codex in this repo and ask
  it to run `make fit`. A reply quoting *"Blocked: that make target is compute"* means the
  hooks load. *"LoginNodeError"* means they do not, and only the in-code guard fired. Either
  way nothing runs. Report the result here.
- **Codex user-level skills overlap.** `~/.codex/skills/` has `conventional-commits` and
  `git-pr-workflow`; the repo ships `conventional-commits` and `git-commit-and-pr`. Inside this
  repo the repo versions are authoritative.
- ~~**`make check` is weak for now.**~~ Superseded 2026-09-29: 171 fast tests and 173 total
  tests pass, including censoring/coverage, model, metric, partition and cache-integrity tests.
  Historical note: tests initially covered only config, schema and guards.
- **PyMC 6.3.2 / arviz 1.3.0** resolved — newer than the reference notebooks (5.28 / 1.1).
  Every API the rules name was verified present; notebook idioms may still not transfer.
- **Final integration checks:** no clean-commit numerical reproduction and no LaTeX-engine
  compilation. Source/artifact integrity, table structure and PDF/PNG layout were checked;
  these do not replace either missing check. Absolute artifact paths require revalidation
  when moving the analysis. See `docs/final-provenance-review.md`.
- **Professor's conformal-method requirement (updated 2026-10-08):** the primary score
  envelope remains an adaptation, and the new post hoc comparator directly implements
  Qin's two-sided Cox bootstrap variant. Candès/Gui have not been implemented. Supervisor
  agreement is still needed before claiming the proposal's method requirement is met;
  Qin's DSS censoring support and molecular refit stability are weak.

## Handoff log

### 2026-10-08 — Codex — Survival-uncertainty reliability goal achieved
- **Changed:** added conservative lower bounds, fixed-ridge Qin refits, training support
  diagnostics and isolated stage-05 orchestration. Added complete paired comparison report,
  four PDF/PNG figures, reproducible export/audit scripts and method/results documentation
  in `docs/uncertainty-reliability.md`; Git/PR drafts updated.
- **Verified:** SLURM 1337 ran `make format`, `make check` (228 passed), `make test-all`
  (232 passed). Real job 1338 completed in 5m57s: all 4,000 ridge refits and 96 interval
  cells. Job 1341 passed final `make check` (228), independent audit (2,523/2,523), and
  aggregate-only exports. Four PNGs visually inspected; `git diff --check` passed.
  Independent lower-score, numerical-method, leakage and statistical reviews found no
  remaining confirmed implementation defect. Original result hashes remain unchanged.
- **Results:** previous Qin refits totaled 3,829/4,000; ridge completes all. At 90%
  nominal/120 months, DSS lower-bound widths shrink by 3.8–8.8 months, but observed
  coverage ranges cross nominal and sets still span 85.5–90.1 months. OS tradeoffs are
  mixed. DSS 300-month IPCW remains unsupported; full Qin event-pool ESS remains 12.8.
  Some PAM50 subgroups retain low realised coverage despite broad intervals.
- **Not done / caveats:** post hoc reused-test comparisons, different one-sided question,
  no individual clinical guarantee, no independent external validation or clean-commit
  numerical reproduction. Jobs 1334–1336 and 1339–1340 stopped on test/formatting issues,
  corrected before the successful runs above. No commits, pushes or PR creation.
- **Next:** use `docs/git-checkpoints.md` for reviewable commits, then review the full
  comparison with the supervisor before considering any new primary method or horizon.

### 2026-10-08 — Codex — Qin bootstrap CPI exploratory comparator
- **Changed:** added Qin et al. Cox bootstrap CPI, isolated stage-05 mode/config, synthetic
  tests and `docs/qin-exploratory.md`. The original evaluation plan and coverage JSON retain
  their original manifest hashes. Job 1329 was canceled during review; its plan is archived.
- **Verified:** independent read-only method/leakage review found no remaining confirmed
  defect after corrections. SLURM 1332: `make check` 176 passed, `make test-all` 179 passed.
  SLURM 1333 completed four arms, 24 interval cells, a source/input-bound exploratory plan
  and aggregate manifest. Original primary plan/coverage SHA-256 values match their manifest.
- **Findings:** 90% DSS 300-month Qin median widths are 272.1/270.3 months for clinical/
  molecular, versus primary 274.4/273.1. DSS 300-month IPCW remains unavailable; 323
  observed DSS events yield only 12.8 inverse-weighted effective event-pool size. The DSS
  molecular bootstrap passed exactly 900/1,000 refits, the preset minimum. See the report.
- **Caveats:** post hoc test reuse, marginal independent-censoring assumption, PH departures,
  weak DSS tail and dirty-tree source archive prevent a confirmatory validity claim.
- **Next:** user commits the existing pipeline checkpoints then this Qin checkpoint, seeks
  supervisor review of exploratory method status, and performs clean-commit reproduction.

### 2026-10-08 — Codex — Aggregate results and model orientation
- **Changed:** documented a read-only sense check for the user's model/metric overview;
  no statistical code, config, patient data or saved results changed.
- **Verified:** reviewed `outputs/reports/thesis_results.md`, aggregate EDA/performance/
  coverage JSON, the frozen evaluation plan, method documents and the prior independent
  audit. SLURM 1325 `make check` passed: Ruff clean, 171 tests passed, 2 slow deselected.
  No new fit or numerical result was computed in this review.
- **Assessment:** results are internally coherent for exploratory internal validation:
  DSS Uno C is about 0.686–0.705; molecular gains are modest; primary 90% conservative
  intervals are about 271–274 months wide at a 300-month horizon. Heavy censoring makes
  exact latent coverage unobservable. Four Cox PH diagnostics flag departures, Weibull
  family/prior sensitivity matters, and DSS tail IPCW is unsupported at 300 months.
- **Caveats:** independent integrity review does not establish clinical utility, external
  validity or clean-commit numerical reproduction. Explain the distinction between
  ranking, probability error, posterior uncertainty and coverage bounds to the user.
- **Next:** user checkpoint commits/PR and thesis integration; obtain supervisor agreement
  on the conservative conformal adaptation before calling the proposal fully matched.

### 2026-09-29 — Codex — Checkpoint 5: thesis exports and final handoff
- **Changed:** completed all six generated RQ answers and the thesis report; exported 137
  PDF/PNG pairs and 94 LaTeX tables; added reusable provenance audit, runbook, review records
  and all five checkpoint command/PR drafts. Closed this session's claim.
- **Verified:** export job 1224 recorded 378 artifacts; independent artifact review 1226
  checked dimensions, numerical consistency and table structure. Independent provenance
  job 1230 passed 906 checks. Integrated SLURM 1231 ran `make format`, `make check`
  (171 passed), `make test-all` (173 passed), and
  `uv run python scripts/audit_results.py --require-exports --output outputs/reports/provenance_audit.json`
  (906 checks, zero failures). Final documentation-only edits followed; no statistical or
  rendering source changed. Independent review confirmed checkpoint staging covers all
  intended code/config/test/docs files and excludes patient artifacts.
- **Review:** independent leakage, conformal, statistical, orchestration and artifact
  findings resolved. Failed initialization and audit-development checks remain documented.
  No model, preprocessing, calibration or statistical evaluation choices changed after test access.
- **Caveats:** this is internal validation with PH departures, broad prediction intervals,
  exploratory subgroup comparisons and unsupported DSS tail IPCW. Integrity checks verify
  archived dirty sources, not clean-commit numerical reproduction. LaTeX compilation and
  optional Codex hook integration remain unverified. Nothing staged, committed or pushed.
- **Next:** user runs `docs/git-checkpoints.md`, integrates exports in the thesis and performs
  clean-commit reproduction; preserve the frozen original analysis and its manifests.

### 2026-09-29 — Codex — Checkpoint 4: frozen final evaluation
- **Verified:** SLURM 1211 completed stage 05 after `make check` (168 passed) and
  `make test-all` (170 passed). Twelve models evaluated on 495 patients, with paired
  1,000-replicate primary-event-stratified bootstrap; every performance replicate succeeded.
  Saved 168 interval-report cells, subgroup summaries and paired uncertainty comparisons.
- **Review:** independent aggregate statistical checks (1217/1219), numerical/table consistency
  review (1218), and provenance audit (1214: 272 checks, zero errors). No invalidating
  statistical finding. Final exports undergo a separate visual/freshness pass.
- **Findings:** modest primary DSS PAM50 discrimination gains for RSF and Weibull; Cox gain
  uncertain. Primary conservative 90% median widths 270.8–274.4 months at a 300-month
  restricted horizon. Favorable marginal realised-coverage bounds do not imply useful
  precision or subgroup guarantees. DSS300 IPCW support is insufficient and stays unavailable.
- **Caveats:** exploratory comparisons, fixed-fit bootstrap, PH departures, net-DSS/competing
  death assumptions and weak tail support are explicit in the generated results report.
- **Next:** finish PDF/PNG/LaTeX and generated RQ report checks, then close the active claim.

### 2026-09-29 — Codex — Checkpoint 3: calibrated restricted survival intervals
- **Changed:** conservative score-envelope CP for restricted event time, stratified finite
  ranks with a worst-stratum threshold, separate approximate IPCW sensitivity, coverage
  bounds/width/subgroup summaries and fixed-primary-event-stratum bootstrap intervals.
- **Verified:** 1210 completed 144 calibration states for 12 fits × three alphas × two
  horizons × two methods; 168 fast tests passed. Independent mathematical review confirmed
  score domination, finite-sample ranks, censoring ties and unknown-stratum prediction;
  bootstrap success-count and stratified-resampling findings fixed and tested.
- **Caveats:** coverage of min(T, horizon), not unrestricted T. IPCW is assumption-dependent;
  unavailable tail support is retained explicitly. Native Bayesian intervals remain separate.
- **Next:** job 1211 passed `make check` (168) and `make test-all` (170), then opened final
  evaluation under a fixed plan. No model/preprocessing/calibration changes after this point.

### 2026-09-29 — Codex — Checkpoint 2: complete training model matrix
- **Changed:** shared model interface, Cox/PH score tests, untuned RSF, Bayesian Weibull
  and lognormal AFT, strict MCMC gate, frozen fit index and training-only sensitivities;
  hypothetical-profile credible bands versus individual predictive intervals.
- **Verified:** SLURM 1204 `make format && make check` (163 passed, 2 slow deselected),
  real data/split/EDA reuse and all 16 fits completed. All eight Bayesian fits passed
  unchanged R-hat/ESS/zero-divergence gates. SLURM 1210 check: 168 passed; eight hypothetical
  profiles generated 16 PDF/PNG files. Calibration continues in the same job.
- **Failure resolved:** 1187 OS molecular chain stuck at initialization (1,000 divergences,
  R-hat 1.529); finite log probability concealed non-finite gradients. Seeded bounded
  prior-centered starts plus finite-gradient preflight fixed it (focused job 1200: R-hat
  1.004, bulk ESS 2,313, tail ESS 2,415, zero divergences). Prior fits archived, then the
  entire matrix refitted with common settings. Priors and diagnostic thresholds unchanged.
- **Review:** independent orchestration review found no leakage; source/protocol binding,
  prediction-cache integrity and gate-status propagation findings fixed. Integration job
  1207 caught source files being added during its cache-recovery test; stable rerun 1210 passes.
- **Caveats:** PH/PPC assumptions must be discussed with results; passing MCMC is not
  evidence of a correct likelihood family. Real test remains unopened; no commits/pushes.
- **Next:** complete calibration, freeze evaluation code and run stage 05 once.

### 2026-09-29 — Codex — Checkpoint 1: real data and training EDA
- **Changed:** validated source loading/cohort preparation; frozen physical partition tables;
  train-only preprocessing and EDA; deterministic source archives; analysis decisions.
- **Verified:** SLURM 1185 ran `make format && make check` (96 passed), stages 01/02 and
  EDA; 1,979 eligible patients, split 989/495/495, 18 figure files. Synthetic split archived.
  SLURM 1186 ran `make format && make check` (149 passed, 2 slow deselected, including
  newly integrated model/conformal unit tests), then reused stage 02 without redrawing.
- **Review:** independent data/leakage review found no endpoint or leakage errors; fixed
  required raw-checksum validation and stage-02 source-archive provenance findings.
- **Caveats:** runs record dirty Git plus immutable source/config archives. Models have not
  yet been fitted to real data; no real test evaluation. Generated files remain ignored.
- **Next:** continue active claim through model fitting, conformal evaluation and exports.

### 2026-09-29 — Codex — Stage 02 frozen stratified splits
- **Changed:** `src/brca/data/splits.py`, `src/brca/manifest.py`, stage-02 CLI, string-ID
  preservation in `load_cohort`, `tests/test_splits.py`, Makefile `SPLIT_ARGS` and corrected
  compute-stage comment, `docs/experimental-design.md` usage/schema/freeze contract.
- **Verified:** `uv sync` → resolved 103 / checked 81 packages; `make format` clean;
  final `make check` → ruff clean, **61 passed**. First focused run had 34 passed / 1 failed:
  `TypeError: Invalid value '123' for dtype 'str'`; fixed the malformed-ID test fixture by
  casting to object before inserting an integer. No product-code failure in that check.
- **Verified:** `make split SPLIT_ARGS='--synthetic --allow-dirty'` succeeded; reran the same
  command through an assertion harness → unchanged split bytes and mtime, fresh manifest
  with `split_reused: true`, matching config/split hashes. Synthetic sizes: 200/100/100;
  event rates: 72.5%/72.0%/72.0% (train/calibration/test). `git check-ignore` confirmed the
  split and manifests are ignored. No real patient data read and no model compute run.
- **Review:** independent static leakage/reproducibility review of an isolated code snapshot;
  fixed its synthetic-symlink labelling finding and added a regression test. No remaining
  findings; reviewer did not execute tests.
- **Not done / caveats:** stage 01 decisions unchanged. Generated `data/processed/splits.json`
  and two `outputs/manifests/02_make_splits_*.json` files are synthetic development artefacts.
  Git has no commits: manifests honestly record null SHA / dirty tree. Nothing staged,
  committed or pushed. Downstream stages remain stubs.
- **Next:** user/supervisor resolve the four stage-01 decisions; implement stage 01 to write
  `cohort.csv`. Archive the synthetic split/manifests (or select separate output directories)
  before drawing the real split. User creates the initial commit containing the scaffold
  and stage 02; commit/push commands were drafted, not executed.

### 2026-09-29 — Claude — Cross-checked data against Simona's thesis
- **Changed:** `docs/metabric-data-dictionary.md` (new § Cross-check; missing counts corrected
  to the analysis cohort instead of all 2,509), Open questions in this file.
- **Verified:** read `references/simona.md` §§ 1–8; re-derived her cohort on our data: same two
  exclusions, 1,979 patients, DSS event rate 32.6% (hers 32.7%), missingness profile matches
  within ~1 point except positive lymph nodes (76 missing here, none reported by her).
- **Not done / caveats:** why Kaggle has 1,904 vs our 1,981 is a hypothesis (patients without
  expression/mutation data), unverified without Kaggle IDs. No code changed.
- **Next:** unchanged — user decides the four stage-01 questions (Simona's choices are now
  listed beside each as precedent).

### 2026-09-29 — Claude — METABRIC clinical data acquired
- **Changed:** `src/brca/data/cbioportal.py`, `scripts/00_download_data.py`, `make download`,
  `data_source` block in `config/default.yaml`, `tests/test_cbioportal.py` (offline),
  `docs/metabric-data-dictionary.md` rewritten from the real data, `AGENTS.md` conventions 2 + 9.
  Raw snapshot in `data/raw/cbioportal_brca_metabric/` (gitignored), checksums in
  `data/raw/CHECKSUMS.sha256`.
- **Verified:** `make download` → 2,509 patients / 2,509 samples; `sha256sum -c` OK; second run
  refuses to overwrite; git sees only `CHECKSUMS.sha256`. Cohort check (aggregates): 1,981 with
  OS + PAM50, 1,815 complete cases, KM subtype ordering as published. `make check` → 24 passed.
- **Not done / caveats:** stage 01 not implemented — blocked on endpoint, PAM50 levels,
  missing-data and admin-censoring decisions (Open questions). Nothing committed.
- **Next:** user answers the four open questions; then implement stage 01 (typing, mapping per
  the data dictionary, exclusions logged as counts). Stage 02 is unblocked meanwhile.

### 2026-09-29 — Claude — Codex compatibility + STATUS.md
- **Changed:** `AGENTS.md` (new, canonical), `CLAUDE.md` (now imports it), `STATUS.md`,
  `status-handoff` skill, `.codex/hooks.json`, `.agents/skills` → `.claude/skills` symlink,
  `.githooks/pre-commit`, `src/brca/compute_guard.py` (wired into stages 03–05),
  `tests/conftest.py`, `tests/test_compute_guard.py`, hook payload parsing for both agents,
  SessionStart hook, Stop-gate STATUS check, `make setup`/`hooks`/`check`,
  `docs/status-archive.md` (empty; for log rotation).
- **Verified:** `make check` → 19 passed, ruff clean. Guards 16/16 on Claude- and Codex-shaped
  payloads. Pre-commit blocks a staged `.csv` (exit 1) and nudges on missing STATUS.md.
  Stage 03 raises `LoginNodeError` here. Slow tests exit 3 on the login node. Fixture
  regenerates when deleted. `core.hooksPath` set to `.githooks`.
- **Not done / caveats:** Codex hook loading unverified (probe under Known issues). Nothing
  committed.
- **Next:** user commits the setup (use the `git-commit-and-pr` skill for the commands).
  Then implement stage 02 `make_splits` against the synthetic fixture, with tests for
  disjointness and determinism.

### 2026-09-20 — Claude — initial agent setup
- **Changed:** `CLAUDE.md`, 12 rules, 10 skills, 4 review agents, 4 hooks, 10 docs, runnable
  skeleton (`pyproject.toml`, `Makefile`, `config/default.yaml`, stubs, synthetic fixture).
- **Verified:** `uv lock`/`uv sync` clean; 13 tests passed; ruff clean; guard and Stop-gate
  cases exercised directly.
- **Next:** superseded by the 2026-09-29 entry.
