# METABRIC data dictionary

> Pulled in by: `onboard-metabric-data` skill; `.claude/rules/data-and-privacy.md`.

## Source and provenance

| | |
| --- | --- |
| Source | cBioPortal REST API, study `brca_metabric` — *Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)* |
| Access | **public study** (`publicStudy: true`); clinical tables only |
| cBioPortal import date | 2026-01-05 |
| Downloaded | 2026-09-29 by `scripts/00_download_data.py` (`make download`) |
| Raw snapshot | `data/raw/cbioportal_brca_metabric/` — gzipped API JSON + `PROVENANCE.json` |
| Checksums | `data/raw/CHECKSUMS.sha256` (committed; `cd data/raw && sha256sum -c CHECKSUMS.sha256`) |
| Cite | Curtis et al. 2012 (PMID 22522925), Pereira et al. 2016 (27161491), Rueda et al. 2019 (30867590), and cBioPortal (Cerami 2012, Gao 2013) |

Expression and raw genomic data remain controlled-access at EGA `EGAS00000000083`; this thesis
does not need them (PAM50 is provided as a label). Public availability does **not** make this
repo a redistribution point: `data/` stays gitignored and patient rows never leave the machine.

`data/raw/` is immutable — `make download` refuses to overwrite an existing snapshot.

## Cohort (aggregates from the 2026-09-29 snapshot)

| | n |
| --- | --- |
| Patients in study (1 sample each, all female) | 2,509 |
| With overall-survival outcome **and** PAM50 | 1,981 |
| Complete cases: all predictors + PAM50 + OS | 1,815 |

The 528 patients without outcome are the 2016 targeted-sequencing extension: no survival, no
PAM50, no PR/HER2. They cannot enter a survival model and are excluded at stage 01. The 1,981
match the published ~1,980-tumour METABRIC cohort.

Outcome summary (n = 1,981): 1,144 deaths (57.7%), of which **646 died of disease and 497 of
other causes**; median follow-up 196 months (reverse Kaplan-Meier); maximum 355 months. Relapse
status is available for 2,488 patients (1,002 recurred).

Kaplan-Meier median overall survival by subtype (months): claudin-low 219, LumA 187, Normal 159,
Basal 131, LumB 123, Her2 107. The ordering of LumA above Basal and Her2 matches the published
prognostic ordering, which is evidence the fields are mapped correctly.

## Source → schema mapping

| Schema column | cBioPortal attribute | Level | Notes |
| --- | --- | --- | --- |
| `patient_id` | `patientId` | patient | unique; 1 sample per patient |
| `survival_months` | `OS_MONTHS` (or `RFS_MONTHS`) | patient | months; depends on the endpoint decision |
| `event` | `OS_STATUS` `1:DECEASED` / `0:LIVING` | patient | → bool; see endpoint below |
| — | `VITAL_STATUS` | patient | Living / Died of Disease / Died of Other Causes |
| `age_at_diagnosis` | `AGE_AT_DIAGNOSIS` | patient | years; 0 missing |
| `tumor_size` | `TUMOR_SIZE` | sample | mm; string in API → float; 25 missing |
| `tumor_grade` | `GRADE` | sample | 1 / 2 / 3; 87 missing |
| `lymph_nodes_positive` | `LYMPH_NODES_EXAMINED_POSITIVE` | patient | string in API → int; 76 missing |
| `er_status` | `ER_STATUS` | sample | Positive / Negative (`ER_IHC` is the IHC-only variant) |
| `pr_status` | `PR_STATUS` | sample | |
| `her2_status` | `HER2_STATUS` | sample | (`HER2_SNP6` is the copy-number variant) |
| `pam50_subtype` | `CLAUDIN_SUBTYPE` | patient | LumA 700, LumB 475, Her2 224, claudin-low 218, Basal 209, Normal 148, NC 6; 1 missing |

Missing counts are within the 1,979-patient analysis cohort (1,981 with OS, minus one missing
cause of death and one with `OS_MONTHS <= 0`). PR, HER2 and PAM50 each have 1 missing. Also
missing but outside the design: `TUMOR_STAGE` 26.0%, `THREEGENE` 11.0%, `LATERALITY` 5.6%,
`CELLULARITY` 3.2%.

Also available but outside the current design: treatment (`CHEMOTHERAPY`, `HORMONE_THERAPY`,
`RADIO_THERAPY`, `BREAST_SURGERY`), `NPI`, `INTCLUST`, `TUMOR_STAGE`, `CELLULARITY`,
`INFERRED_MENOPAUSAL_STATE`, `COHORT` (recruitment site 1–9). Treatment variables are
post-diagnosis decisions; including them as predictors needs thought (they partly encode the
clinician's own prognosis).

## Cross-check against Simona's thesis (`references/simona.md`)

Her thesis (Macedonian; Statistical Modelling course) analysed the same METABRIC study, from a
**different distribution**: the Kaggle re-upload (`raghadalharbi/breast-cancer-gene-expression-
profiles-metabric`), 1,904 patients with 489 expression z-scores and 173 mutation columns. Ours
is cBioPortal's official public study, clinical tables only, 1,981 patients with outcome.

Reproducing her cohort definition on our data (checked 2026-09-29):

| | Simona (Kaggle) | Ours (cBioPortal) |
| --- | --- | --- |
| Exclusions | 1 missing cause of death, 1 time = 0 | exactly the same two |
| Cohort after exclusions | 1,902 | 1,979 |
| Breast-cancer-specific event rate | 32.7% | 32.6% (646 events) |
| Missing `TUMOR_STAGE` / `THREEGENE` / laterality | 26.3% / 10.7% / 5.5% | 26.0% / 11.0% / 5.6% |
| Missing grade / cellularity / tumour size | 3.7% / 2.8% / 1.0% | 4.4% / 3.2% / 1.3% |
| Missing positive lymph nodes | none reported | **76 (3.8%)** |

The ~77 extra patients here are most likely those without expression/mutation profiles, which
Kaggle omits; that would also explain the lymph-node gap. **Unverified** — confirming it needs
the Kaggle patient IDs (both use `MB-xxxx`).

Her choices, as precedent for our open decisions:

- **Endpoint: breast-cancer-specific survival** — death from other causes censored. Overall
  survival re-run only as a sensitivity analysis; Fine-Gray competing risks exploratory only.
- **PAM50: "PAM50 + Claudin-low"**, claudin-low kept as its own level. NC is not discussed.
- **Missing data:** `TUMOR_STAGE` dropped (26% missing; its information is carried by size,
  nodes and grade). Tumour size by MICE; grade and cellularity by mode, after MICE collapsed
  their category distributions. She names an explicit `unknown` level as the better option for
  ordinal variables. Her imputation section precedes the train/test split, and only gene
  selection and scaling are stated as train-only, so imputation **appears** to have been fit on
  the full data — our `splits-and-leakage` rule requires train-only.
- **Treatment variables excluded** (pre-treatment prognosis; confounding by indication).
- **No administrative truncation**; IBS computed only up to **120 months**, because Brier is
  unstable once few patients remain at risk.
- **NPI as a clinical benchmark** (C-index 0.6796 on her test set). NPI is in our data too.
- **Split 75:25**, stratified by event. Ours is 50/25/25 because conformal needs a calibration set.

Two of her findings bear directly on our models:

- **Weibull fitted worst of the three AFT families** on the clinical-molecular features (AIC:
  log-normal 6272.7 < log-logistic 6307.5 < Weibull 6342.0; C 0.7033 / 0.7023 / 0.6965). The
  brief asks for a Bayesian Weibull/AFT, so a log-normal likelihood deserves at least a
  posterior-predictive comparison — a misfitting likelihood gives miscalibrated Bayesian
  intervals, which is exactly what RQ5–RQ6 probe.
- **ER status violates proportional hazards** in Cox (Schoenfeld statistic 15.87, p < 0.0005);
  she used ER-stratified Cox as a robustness check.

Her "clinical-molecular" group already **includes** PAM50 and her contrast was
clinical-molecular vs +genomics, so her results do not answer our RQ2. Her unstratified Cox on
clinical-molecular features (C 0.7005) is a useful reference point for ours.

## Decisions still open — they change stage 01

1. **Endpoint.** Overall survival treats all 1,144 deaths as events. Disease-specific survival
   treats the 497 other-cause deaths as censored, which assumes they are uninformative — a
   competing-risk assumption to state explicitly. Relapse-free survival uses `RFS_*`.
   *Simona used disease-specific; following her keeps the two theses comparable.*
2. **PAM50 levels.** The schema expects 5 levels. Claudin-low (218) must be its own level or
   folded; NC (6) is too small to model and is best mapped to `unknown` or dropped, identically
   in both feature-set arms. *Simona kept claudin-low as its own level.*
3. **One patient with `OS_MONTHS <= 0`.** The schema rejects it by design; exclude it at stage 01
   with an explicit logged count.
4. **Missing predictors** (166 of the 1,981 lack at least one; mostly grade 87, nodes 76,
   size 25): complete-case analysis versus imputation fit on train only. Either way, identical
   in both arms. *Simona imputed (MICE / mode); we must fit any imputer on train only.*
5. **Administrative censoring.** `config.outcome.administrative_censoring_months` is 300, but
   follow-up runs to 355. *Simona did not truncate and capped IBS at 120 months instead.*
