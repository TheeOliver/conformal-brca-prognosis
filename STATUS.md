# STATUS

Shared progress ledger for **Claude Code and Codex**. Read it at the start of every session;
update it at the end of every session that changed files. Procedure:
`.claude/skills/status-handoff/SKILL.md`. Newest handoff entries go at the **top** of the log.

## Snapshot

- **Phase:** 1 — stage 02 complete on synthetic data; stage 01 blocked on four decisions
- **Last updated:** 2026-09-29 by Codex
- **Next up:** user decides endpoint, PAM50 levels, missing-data handling and admin censoring;
  then implement stage 01 against `data/raw/cbioportal_brca_metabric/`. Archive the synthetic
  split/manifests before freezing a real-data split. Setup and stage 02 are still uncommitted.

## Active claims

Claim before starting non-trivial work; remove your row when you hand off.

| Task | Agent | Branch | Since |
| --- | --- | --- | --- |
| — | — | — | — |

## Milestones

**Phase 0 — Setup**
- [x] Agent configuration: `AGENTS.md` (canonical), `CLAUDE.md` (imports it), rules, skills,
      review checklists, hooks for Claude and Codex
- [x] Agent-independent guards: `.githooks/pre-commit`, `src/brca/compute_guard.py`
- [x] Environment: uv, Python 3.12, `uv.lock`
- [x] Synthetic fixture + scaffold tests (`make check` → 19 passed)
- [ ] First commit — *user*

**Phase 1 — Data**
- [x] METABRIC clinical data acquired: cBioPortal public study, 1,981 usable patients
      (`make download`; `docs/metabric-data-dictionary.md`)
- [ ] Endpoint decided: overall vs disease-specific survival vs relapse (see Open questions)
- [ ] Stage 01 `prepare_data` + tests
- [x] Stage 02 `make_splits` + tests (synthetic frozen `splits.json`, disjointness asserted)

**Phase 2 — Models** (each for *both* feature sets)
- [ ] Cox PH, with proportional-hazards check
- [ ] Bayesian Weibull AFT: prior predictive, diagnostics gate, PPC, prior sensitivity
- [ ] Random Survival Forest (untuned, by design)

**Phase 3 — Conformal**
- [ ] Censoring adaptation chosen and justified (see Open questions)
- [ ] Split CP with finite-sample quantile + coverage unit tests
- [ ] `conformal-validity-auditor` review — by the agent that did *not* write it

**Phase 4 — Evaluation**
- [ ] Uno C, Brier, IBS, td-AUC with bootstrap CIs on a shared interior time grid
- [ ] Coverage + width at 80/90/95%, marginal and by PAM50 / ER subgroup

**Phase 5 — Answers** (each links to its `outputs/metrics/*.json`)
- [ ] RQ1 · [ ] RQ2 · [ ] RQ3 · [ ] RQ4 · [ ] RQ5 · [ ] RQ6

**Phase 6 — Thesis artefacts**
- [ ] Figures and tables exported per `docs/figure-style.md`
- [ ] `leakage-auditor`, `stats-reviewer`, `reproducibility-checker` all clean

## Blockers & open questions

- ~~**BLOCKER — METABRIC access.**~~ Resolved 2026-09-29: cBioPortal `brca_metabric` is a
  public study with all needed clinical fields + PAM50 + OS/RFS. *(Claude)* Simona used a
  different distribution of the same study (Kaggle, 1,904 patients); her cohort definition
  reproduces on ours (same 2 exclusions, 32.6% vs 32.7% DSS event rate). See
  `docs/metabric-data-dictionary.md` § Cross-check.
- **Endpoint — blocks stage 01.** 497 of 1,144 deaths (43%) are from other causes, so overall
  vs disease-specific survival will give materially different models; relapse-free is also
  available. Decide with the supervisor. **Precedent: Simona used disease-specific** (OS only
  as sensitivity). *Owner: user.*
- **PAM50 levels — blocks stage 01.** Data has claudin-low (218) and NC (< 10) beyond the 5
  schema levels. Proposal: claudin-low as its own level (**as Simona did**), NC (6) →
  `unknown`. *Owner: user.*
- **Missing predictors — blocks stage 01.** 166 of 1,981 lack a predictor (grade 87, nodes 76,
  size 25). Complete-case or train-only imputation? Simona imputed (MICE for size, mode for
  grade), apparently before her split — ours must be train-only. *Owner: user.*
- **Admin censoring.** Config says 300 months; follow-up runs to 355. Simona did not truncate
  and capped IBS at 120 months. Confirm. *Owner: user.*
- **AFT likelihood (Phase 2).** Simona found Weibull the worst-fitting AFT family by AIC
  (log-normal best). The brief asks for Bayesian Weibull/AFT — compare against log-normal by
  posterior predictive checks before committing. *Owner: user + whoever builds the model.*
- **Conformal censoring adaptation.** Candès-style lower predictive bound, IPCW-weighted scores,
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

## Known issues

- **Codex hooks unverified.** `.codex/hooks.json` targets Codex 0.153.4's hook format, but
  could not be exercised from this session. 30-second probe: open Codex in this repo and ask
  it to run `make fit`. A reply quoting *"Blocked: that make target is compute"* means the
  hooks load. *"LoginNodeError"* means they do not, and only the in-code guard fired. Either
  way nothing runs. Report the result here.
- **Codex user-level skills overlap.** `~/.codex/skills/` has `conventional-commits` and
  `git-pr-workflow`; the repo ships `conventional-commits` and `git-commit-and-pr`. Inside this
  repo the repo versions are authoritative.
- **`make check` is weak for now.** Tests cover config, schema and guards — not the science.
  It becomes meaningful when conformal coverage tests land (Phase 3).
- **PyMC 6.3.2 / arviz 1.3.0** resolved — newer than the reference notebooks (5.28 / 1.1).
  Every API the rules name was verified present; notebook idioms may still not transfer.

## Handoff log

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
