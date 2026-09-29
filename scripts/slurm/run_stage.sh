#!/bin/bash
# Submit one pipeline stage to SLURM.
#
#   sbatch scripts/slurm/run_stage.sh 03
#
# This exists because stages 03-05 must never run on the login node -- see
# CLAUDE.md convention 1 and docs/slurm-guide.md. The PreToolUse hook that
# refuses those commands points here.
#SBATCH --job-name=brca-stage
#SBATCH --partition=openlab-queue
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm-%j.out

set -euo pipefail

STAGE="${1:?usage: sbatch scripts/slurm/run_stage.sh <stage NN, e.g. 03>}"
CONFIG="${2:-config/default.yaml}"

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

# 4 cores, one per MCMC chain (config.mcmc.chains). Keep BLAS single-threaded so
# the chains do not oversubscribe the allocation fighting each other for cores.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

shopt -s nullglob
matches=(scripts/"${STAGE}"_*.py)
if [ ${#matches[@]} -ne 1 ]; then
  echo "expected exactly one script for stage '${STAGE}', found: ${matches[*]:-none}" >&2
  exit 2
fi

echo "job ${SLURM_JOB_ID:-local} :: stage ${STAGE} :: ${matches[0]} :: $(date -Is)"
uv run "${matches[0]}" --config "${CONFIG}"
