# Reproducing the thesis analysis

All tests and pipeline execution in this session run through SLURM. Prepare the environment
with `make setup`, download the public clinical tables with `make download`, and create the
logs directory before submission. Never stage generated patient or posterior files.

For a new analysis directory, with the registered configuration and a committed source tree:

```bash
mkdir -p logs
sbatch scripts/slurm/run_pipeline.sh config/default.yaml
```

The workflow runs checks, prepares data, freezes the split, explores training data, fits the
12 primary endpoint/model/feature combinations and four DSS Bayesian sensitivity fits,
calibrates, evaluates and exports. It deliberately stops if a Bayesian diagnostic gate fails.
The synthetic split from initial development was archived before the first real split;
existing real splits are reused exactly. Never delete or replace one to seek better results.

To execute a specific stage:

```bash
sbatch scripts/slurm/run_stage.sh check
sbatch scripts/slurm/run_stage.sh 03
sbatch scripts/slurm/run_stage.sh 04
sbatch scripts/slurm/run_stage.sh 05
sbatch scripts/slurm/run_stage.sh 06
sbatch scripts/slurm/run_stage.sh test-all
sbatch --partition=openlab-queue --cpus-per-task=1 --mem=2G --time=00:05:00 \
  --output=logs/audit-%j.out --wrap='uv run python scripts/audit_results.py --require-exports'
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

Explicit development runs can append `config/default.yaml --allow-dirty`. Each records
Git state plus a compressed source/config archive addressed by its SHA256. Checkpoint
commits remain the user's responsibility. Match the manifest's archive/hash, configuration,
raw checksums and `uv.lock` when reproducing numbers; a later Git commit alone is insufficient.

Stage03 reuses only artifacts whose hashes, training data, configuration and fitting source
match. Stage04 uses the saved train-only preprocessing and models. Stage05 writes its frozen
plan before opening test, then stores an ignored integrity-checked test cache and predictions.
After this point stages03/04 refuse to refit/recalibrate. An interrupted stage05 can resume
with the same code and plan; it cannot change the analysis after observing test outcomes.
Stage06 regenerates figures from aggregate JSON and does not open patient tables.

If a fresh experiment is scientifically justified, use a distinct configured output directory
and document the reason before inspecting further test results. This is not permission to
select hyperparameters or methods using the held-out results from the first analysis.

Primary outputs:

- `outputs/metrics/training_eda.json`: predictor, missingness and survival descriptions.
- `outputs/models/fit_index.json`: frozen fits, preprocessing, training time grids and hashes.
- `outputs/metrics/bayes_diagnostics_*.json`: MCMC gates, prior/PPC and coefficients.
- `outputs/metrics/bayesian_sensitivity.json`: training-only alternative-prior/family comparisons.
- `outputs/models/conformal_calibration.json`: calibration ranks and support diagnostics.
- `outputs/metrics/performance.json`: Uno C, IBS, Brier and dynamic AUC with paired bootstrap.
- `outputs/metrics/conformal_coverage.json`: observable coverage bounds, IPCW sensitivity,
  widths and subgroup summaries at all nominal levels and both restricted horizons.
- `outputs/metrics/rq2_feature_set_comparison.json`: paired molecular-minus-clinical differences.
- `outputs/figures/`, `outputs/tables/`: PDF/PNG figures and booktabs LaTeX tables.

IPCW coverage can be unavailable at the 300-month horizon if training censoring survival
falls below the registered support threshold. Report that limitation together with
observable coverage bounds and width. Do not switch the primary horizon in response.

LaTeX exports use `booktabs` and `siunitx`; long coverage appendices also use `longtable`.
Include these packages in the thesis preamble. The cluster has no `pdflatex`, so rendering
checks verify generated content and figures; compile the tables in the thesis document.
