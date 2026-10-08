# Final orchestration and provenance review

Review date: 2026-09-29. The review was independent of the root agent's pipeline
orchestration. Conformal mathematics received a separate review by the classical
model agent; this note covers caller integration and saved-result integrity.

**Outcome: artifact integrity checks passed, with qualified reproducibility.**
SLURM job **1230** completed **906 checks with no failures** after export job 1224.
The machine-readable record is `outputs/reports/provenance_audit.json`; its audit
script digest is recorded inside the report. No patient tables, prediction caches,
or posterior/model objects were decoded by this audit.

## Verified integration

- Models and preprocessing fit training patients only. Stage04 loads training and
  calibration partitions. Stage05 opens the test partition after publishing its
  frozen plan and binds cached predictions to that plan and model artifacts.
- Both DSS and OS use the primary DSS event labels for calibration strata and
  bootstrap strata. Clinical and clinical+PAM50 comparisons share patient draws.
- Native Weibull predictive intervals remain separate from conformal intervals.
  Restricted survival uses the registered 300- and 120-month horizons; full model
  follow-up is retained. Unsupported IPCW calibration produces full-support
  intervals with an explicit support status, not an asserted coverage guarantee.
- The current calibration file matches the successful stage04 manifest and the
  frozen stage05 plan. Stage06 requires that same producing calibration manifest.
  Frozen analysis-source and prediction-cache bindings prevent silent analysis
  changes on resumed evaluation.

## Executed checks

| Producing stage | SLURM job | Verified artifacts |
| --- | --- | --- |
| Stage03 fitting | 1204 | Fit index plus 43 registered artifact digests; 12 main and 4 sensitivity fits |
| Stage04 calibration | 1210 | Calibration digest; all 144 registered model/alpha/horizon/method combinations |
| Stage05 evaluation | 1211 | 16 manifest artifact digests and frozen plan bindings |
| Stage06 export | 1224 | 378 output digests, 22 input digests, and 7 current rendering-source files |

Three required source archives matched their SHA-256 digests. Configuration and
lockfile bytes matched the archived copies; each producing manifest recorded 83
package versions consistent with `uv.lock`. All 27 frozen analysis-file digests
matched. Metadata links the processed cohort to a `real_metabric` preparation
manifest and the frozen split; patient-table bytes were not reread for this check.

Calibration diagnostics contain 72 primary score-envelope states, 54 available
IPCW sensitivity states and 18 states with unavailable censoring support. Those
status counts describe calibration availability, not measured latent coverage.

The reusable audit's synthetic tests detect altered or missing artifacts, mismatched
calibration/config digests, and source archives containing disallowed data paths or
path traversal. `make check` in the isolated audit worktree passed in job 1230:
Ruff clean, **106 tests passed**. This is that worktree's test count; the integrated
pipeline's test count is recorded separately by the main agent.

Earlier audit-development jobs exposed summary/schema assumptions and lint issues:
1213 used the wrong calibration status location; 1220/1225 failed lint; 1229 expected
an optional stage01 status field. These were corrected before job 1230. None was a
failure of a fitted artifact or an observed patient result.

## Reproducibility qualifications

The producing jobs explicitly record `git_dirty: true`, Git baseline
`d22268b40aa9cd86411a0a0519705fef918ab557`, and complete source/config/lock archives.
This verifies the archived code bytes; it does **not** satisfy the repository's
strict clean-commit reproduction requirement. Committing remains the user's task.
No second numerical fit/evaluation was performed by this audit, so numerical
reproduction and byte-identical MCMC reruns are not claimed.

Source archives include runnable `src/`, `scripts/`, `config/`, `pyproject.toml` and
`uv.lock`; they do not substitute for committing the complete test/documentation
tree. Manifests and fit indices contain absolute paths: restoration elsewhere needs
the original directory layout or an explicit, revalidated path migration.

To repeat the integrity check after copying the audit script into the main tree:

```bash
sbatch --job-name=brca-audit --partition=openlab-queue --cpus-per-task=1 \
  --mem=2G --time=00:05:00 --output=logs/audit-%j.out \
  --wrap='uv run --frozen python scripts/audit_results.py --require-exports'
```

The audit requires SLURM and exits nonzero if an integrity check fails. Its report is
written separately from the frozen statistical outputs.
