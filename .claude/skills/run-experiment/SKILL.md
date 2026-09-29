---
name: run-experiment
description: Run pipeline stages (prepare data, split, fit models, calibrate conformal, evaluate) on the SLURM cluster with a reproducible manifest. Use whenever asked to run the pipeline, fit or refit the models, rerun a stage, regenerate results, produce metrics, or reproduce a previous run.
---

# Run an experiment

The pipeline is five stages. Stages 03-05 are **compute** and must not run on the login node.

## 1. Establish scope

Ask which stages are needed, or infer from what changed:

| Changed | Rerun from |
| --- | --- |
| raw data, schema, preprocessing | `01_prepare_data` |
| seed or split fractions | `02_make_splits` (invalidates **everything** downstream) |
| a model, priors, RSF settings | `03_fit_models` |
| nonconformity score, alphas | `04_conformal` |
| metrics, time grid, subgroups | `05_evaluate` |

Never rerun `02` casually — it redraws the split and silently invalidates every cached fit
and every calibrated interval. Confirm with the user before doing so.

## 2. Pre-flight

1. `git status` — a dirty tree means the run is untraceable. Commit or state clearly that
   this is a throwaway run that will not be reported.
2. Confirm `config/default.yaml` holds the intended seed, split and alphas. Read it; do not
   assume.
3. `uv sync` so the environment matches `uv.lock`.
4. `make test` — if the fast tests fail, the pipeline output is not worth producing.

## 3. Choose the execution mode

Read `docs/slurm-guide.md` and pick:

- **Stages 01-02** are cheap; run them directly with `uv run python scripts/0N_*.py`.
- **Stages 03-05** go through `sbatch`. Sampling four MCMC chains needs four cores; ask for
  them explicitly rather than accepting the default. A PreToolUse hook will refuse these on
  the login node — that refusal is correct, do not work around it.

## 4. Submit and monitor

```bash
sbatch scripts/slurm/run_stage.sh 03           # returns a job id
squeue -u "$USER"                              # queued / running
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS   # after it ends
```

Report the job id back to the user immediately, then poll rather than blocking.

## 5. Verify before believing

1. Job state is `COMPLETED`, not `FAILED` or `TIMEOUT`. Read the SLURM log even on success —
   a stage can exit 0 having skipped work.
2. A manifest exists in `outputs/manifests/` with the expected git SHA and config hash.
3. For stage 03, the MCMC diagnostics gate passed (`.claude/rules/mcmc-diagnostics.md`). A
   fit that failed the gate is not a result.
4. Expected files are present in `outputs/metrics/`, and non-empty.

## 6. Report

State which stages ran, the job ids, where outputs landed, and anything that failed the gate.
If the run used `tests/fixtures/synthetic_metabric.csv` rather than real METABRIC, **say so
prominently** — synthetic numbers are never presented as findings.

Design and rationale: `docs/experimental-design.md`. Cluster details: `docs/slurm-guide.md`.
