# Running on the cluster

> Pulled in by: `run-experiment` skill; referenced by the login-node PreToolUse hook and
> CLAUDE.md convention 1.

## The rule

`openlab-slurm.hpc.local` is a **login node**: 8 cores, 15 GB, shared by the whole group. It is
for editing, linting, fast tests and job submission. It is not for sampling four MCMC chains or
fitting a 500-tree forest — doing so degrades the machine for everyone.

Stages 03-05, `make all`, `make test-all` and any PyMC sampling go through `sbatch` or `srun`.
A PreToolUse hook enforces this. If it refuses a command, the refusal is correct: submit the
job instead of working around the guard.

## Partitions

| Partition | Time limit | Notes |
| --- | --- | --- |
| `openlab-queue` | 1 day | default |
| `priority-queue` | 7 days | for long runs |
| `cybermacs-queue` | 4 hours | short jobs |
| `vezilka-queue` | 7 days | |

Check availability with `sinfo`.

## Sizing

- **MCMC**: one core per chain. `config.mcmc.chains` is 4, so request `--cpus-per-task=4` and
  set the sampler's `cores` to match. Requesting fewer serialises the chains; requesting more
  wastes an allocation.
- **RSF**: `n_jobs` should equal `--cpus-per-task`, not `-1` — `-1` grabs every core on the
  node, including ones not allocated to the job.
- **Memory**: METABRIC is small (~2000 patients). 8 GB is generous. InferenceData for 4 chains
  x 1000 draws is modest; if memory is the binding constraint, something is wrong.
- **No GPU needed.** Neither PyMC NUTS on this model size nor scikit-survival benefits.

## Submission template

`scripts/slurm/run_stage.sh` takes the stage number:

```bash
sbatch scripts/slurm/run_stage.sh 03
```

It pins BLAS to one thread per chain so the four MCMC chains do not fight for cores:

```bash
#!/bin/bash
#SBATCH --job-name=brca-stage
#SBATCH --partition=openlab-queue
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-%j.out

cd "$SLURM_SUBMIT_DIR"
uv run python "scripts/0${1}"_*.py --config config/default.yaml
```

Record the `SLURM_JOB_ID` in the run manifest so a result can be traced back to its job.

## Monitoring

```bash
squeue -u "$USER"                                            # queued / running
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS         # after it ends
scontrol show job <jobid>                                    # why is it pending
```

Read the log even on success — a stage can exit 0 having skipped its work.
