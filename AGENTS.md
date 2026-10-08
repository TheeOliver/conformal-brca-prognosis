# conformal-brca-prognosis — instructions for every agent

**Claude Code and Codex both work in this repo.** This file is the single source of truth for
project conventions: Codex reads it natively and `CLAUDE.md` imports it. Change conventions here.

Prognostic survival modelling on METABRIC — Cox PH, Bayesian Weibull AFT, Random Survival
Forest — wrapped in split conformal prediction to test whether declared coverage is achieved.

## Session protocol — `STATUS.md` is the shared memory

The two agents alternate and neither remembers the other's sessions. `STATUS.md` is how each
one learns where things stand. Skipping it means redoing or undoing the other agent's work.

1. **Start** — read `STATUS.md` before anything else. Check *Active claims*: never start work
   another agent has claimed. Confirm you are on the branch that work belongs to.
2. **Claim** — before non-trivial work, add a row to *Active claims* (task, agent, branch, date).
3. **Finish** — tick milestones, log decisions, add a *Handoff log* entry at the top saying what
   changed, what you verified (commands actually run, with results), and what is next. Then
   remove your claim.
4. Never delete or rewrite another agent's entries. Append, or mark an entry superseded.

Full procedure: `.claude/skills/status-handoff/SKILL.md`. One agent per working tree at a time;
for parallel work use separate branches or `git worktree`.

## Commands

| Command | What |
| --- | --- |
| `make setup` | `uv sync` + install git hooks (`.githooks/`) |
| `make check` | **definition of done**: ruff + fast tests |
| `make test` / `make test-all` | fast tests / all incl. `slow` MCMC tests (`sbatch` only) |
| `make format` | ruff autofix — run it after editing Python |
| `make download` | stage 00: fetch public METABRIC clinical tables (network only) |
| `make data split fit conformal eval` | pipeline stages 01→05 (03–05 are **compute**) |

Always `uv run <cmd>`, never bare `python`/`pytest` — system Python 3.10 has no dependencies.

## Layout

| Path | Contents |
| --- | --- |
| `src/brca/{data,models,conformal,evaluation,viz}/` | all importable logic |
| `scripts/01…05_*.py`, `scripts/slurm/` | pipeline stages (thin CLI wrappers); SLURM submit |
| `config/default.yaml` | every seed, split fraction, α and threshold |
| `tests/`, `tests/fixtures/` | tests; synthetic fixture generated on first run |
| `data/`, `outputs/` | gitignored; generated, never authored |
| `docs/` | reference material; start at `docs/experimental-design.md` |
| `references/literature/` | 19 full-text papers (up to 3k lines) — grep, don't read whole |
| `.claude/rules/`, `.claude/skills/`, `.claude/agents/` | rules, procedures, review checklists — **for both agents** |

## Hard conventions

1. **Never run compute on this host.** `openlab-slurm.hpc.local` is a shared SLURM login node.
   Stages 03–05, `make all`, `make test-all`, any PyMC sampling → `sbatch`/`srun`
   (`docs/slurm-guide.md`). Stages 03–05 raise `LoginNodeError` here by design.
2. **Patient data never enters git.** The cBioPortal clinical tables are public, but this repo
   is not a redistribution point (expression data stays controlled-access at EGA). No `data/`, `outputs/`, `.csv`/`.parquet`/`.nc`, or notebook outputs in any commit — the
   git pre-commit hook blocks them. Never paste patient-level rows anywhere.
3. **One split, drawn once** by stage 02 into `data/processed/splits.json`. Fit on train;
   conformal calibration on calibration; **test is read only in `scripts/05_evaluate.py`**.
4. **Both feature sets use identical splits.** `clinical` vs `clinical_molecular` is RQ2.
5. **`alpha` is miscoverage** (`0.10` → 90%). Coverage is never reported without width.
6. **Prognostic, not causal.** Never "effect of" / "causes" / "reduces risk" about a
   covariate — write "associated with" / "predictive of". Code, docs and figure labels alike.
7. **`references/notebooks_reference/*.md` are lossy** conversions (truncated dicts, dropped
   arguments). Read for approach; never copy code from them.
8. **The conformal layer is hand-written** (no Python library handles censoring). Every
   nonconformity score needs a unit test asserting coverage.
9. **The pipeline uses the real clinical cohort.** `make download` fetched the
   public cBioPortal clinical tables into `data/raw/cbioportal_brca_metabric/` (immutable;
   `docs/metabric-data-dictionary.md`). Stage 01 prepares the registered DSS/OS cohort;
   the development split is archived. Tests use synthetic fixtures, never real patients,
   and synthetic numbers are never results.
10. **Committing and pushing are the user's.** Draft commands and messages
    (`git-commit-and-pr` skill); never run `git commit`, `git push` or `gh pr create` unless
    told to in that turn. Commit messages and PR text never carry AI attribution.

## Read before editing (Claude loads these automatically; Codex must open them)

| When editing | Read first (in `.claude/rules/`) |
| --- | --- |
| `src/brca/data/**`, `scripts/**`, `data/**` | `data-and-privacy`, `splits-and-leakage`, `reproducibility` |
| `src/brca/models/**` | `sksurv-api`, `survival-labels`; Bayesian: + `bayesian-pymc`, `mcmc-diagnostics` |
| `src/brca/conformal/**` | `conformal`, `splits-and-leakage` |
| `src/brca/evaluation/**` | `evaluation-metrics`, `sksurv-api`, `scientific-claims` |
| `src/brca/viz/**`, `outputs/**` | `figures-and-tables`, `scientific-claims` |
| `notebooks/**` | `notebooks` |
| `docs/**` | `scientific-claims` |
| `tests/**` | `survival-labels` |

## Procedures (`.claude/skills/<name>/SKILL.md`; Codex also finds them via `.agents/skills`)

`status-handoff` (every session) · `run-experiment` · `add-survival-model` ·
`bayesian-model-check` · `conformal-coverage-report` · `compare-feature-sets` ·
`export-thesis-figure` · `onboard-metabric-data` · `distill-paper-note` ·
`conventional-commits` · `git-commit-and-pr`

## Reviews (`.claude/agents/<name>.md`)

`leakage-auditor` · `conformal-validity-auditor` · `stats-reviewer` · `reproducibility-checker`.
Run the relevant one before any result is reported. Claude runs them as subagents; Codex
follows the file as a review checklist, or hands it to a subagent. The reviewer should not be
the agent that wrote the code under review — alternating agents is an advantage here; use it.

## Enforcement — what holds regardless of agent

| Guard | Covers | Where |
| --- | --- | --- |
| data never committed | every committer | `.githooks/pre-commit` (`make setup` installs it) |
| no compute on login node | every runner | `src/brca/compute_guard.py`, `tests/conftest.py` |
| early warnings, auto-format, Stop gate, STATUS snapshot | Claude (verified); Codex (best-effort) | `.claude/settings.json`, `.codex/hooks.json` |

## Definition of done

`make check` passes **and** `STATUS.md` has a handoff entry for the work. Report failures with
their output — never describe unrun checks as passing.
