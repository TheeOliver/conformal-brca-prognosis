---
name: onboard-metabric-data
description: Bring a real METABRIC extract into the project - place it, validate the schema, checksum it, and switch the pipeline off synthetic data. Use when the METABRIC data arrives or access is granted, when asked to load, import or point the project at the real dataset, or when a schema or column-name mismatch appears on load.
---

# Onboard the real METABRIC data

The pipeline currently runs on `tests/fixtures/synthetic_metabric.csv`. This is the one-time
switch to real data. `.claude/rules/data-and-privacy.md` applies to every step.

## 0. Confirm the access position first

METABRIC is controlled-access: EGA `EGAS00000000083`, Synapse `syn1688369`. Before touching
anything, confirm with the user which route granted access and whether any redistribution or
retention terms constrain where the data may live. Do **not** attempt to download it
autonomously, and never fetch it from an unofficial mirror.

`docs/metabric-data-dictionary.md` records the provenance and the expected variables.

## 1. Place it

1. Copy into `data/raw/`, which is gitignored and immutable. Never edit a file there.
2. `sha256sum` every file into `data/raw/CHECKSUMS.sha256` — the one committed file under
   `data/`. This is what makes "which snapshot produced this number" answerable in a year.
3. Record in `docs/metabric-data-dictionary.md`: source, access date, version, file list.

## 2. Validate against the schema, do not adapt to the data

Run the loader. It will fail on mismatch — that is the point.

For each mismatch, decide deliberately:

- **A column is named differently** → add an explicit rename map in `src/brca/data/load.py`.
  Never rename the raw file.
- **A column is missing** → does the analysis change? Adjust `config.feature_sets` and say so
  out loud; silently dropping a predictor changes the research question.
- **A level is unexpected** (e.g. a PAM50 `NC` / `Claudin-low` category) → decide and document
  whether it is its own level or `unknown`. Do not collapse it quietly.
- **Types disagree** (event as `1`/`0` or `"DECEASED"`) → convert in the loader to bool;
  see `.claude/rules/survival-labels.md`.

## 3. Sanity-check the cohort before modelling

Compare against the published METABRIC description — roughly 1,980 tumours:

- n, event rate, median follow-up, censoring proportion;
- a Kaplan-Meier curve overall and by PAM50 — basal-like should track worse than luminal A.
  If it does not, something is mismapped;
- missingness per column, and whether it is plausibly at random;
- `time <= 0`, duplicate patient IDs, impossible ages.

If the cohort does not resemble the published description, stop and report. Modelling a
mis-parsed dataset wastes far more time than checking does.

## 4. Switch over and rerun

1. Point `config.paths.raw` at the real extract.
2. `make data split` — this draws the real split for the first time; it is now frozen.
3. Rerun stages 03-05 via `sbatch` (`run-experiment` skill). Every prior result was synthetic
   and must be regenerated.
4. Run the `leakage-auditor` agent — real data has missingness patterns the synthetic fixture
   does not, and imputation is where leakage enters.

## 5. Keep the synthetic fixture

Do not delete or repoint it. Tests must keep running without real data, both so CI works and
so no test can ever read a patient record.
