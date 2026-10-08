#!/bin/bash
#SBATCH --job-name=brca-pipeline
#SBATCH --partition=openlab-queue
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pipeline-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BRCA_CONFIG="${1:-config/default.yaml}"
if [ "$#" -gt 0 ]; then shift; fi
uv sync --frozen
make check
uv run python scripts/01_prepare_data.py --config "$BRCA_CONFIG" "$@"
uv run python scripts/02_make_splits.py --config "$BRCA_CONFIG" "$@"
uv run python scripts/eda.py --config "$BRCA_CONFIG" "$@"
uv run python scripts/03_fit_models.py --config "$BRCA_CONFIG" "$@"
uv run python scripts/training_uncertainty.py --config "$BRCA_CONFIG" "$@"
uv run python scripts/04_conformal.py --config "$BRCA_CONFIG" "$@"
uv run python scripts/05_evaluate.py --config "$BRCA_CONFIG" "$@"
uv run python scripts/06_export_results.py --config "$BRCA_CONFIG" "$@"
