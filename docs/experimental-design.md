# Experimental design

> Pulled in by: `run-experiment`, `compare-feature-sets` skills;
> `.claude/rules/splits-and-leakage.md`.

## The three-way split

`scripts/02_make_splits.py` draws it **once** from `config.seed`, stratified by event
indicator, and freezes patient IDs into `data/processed/splits.json`.

| Set | Fraction | Used for | Never used for |
| --- | --- | --- | --- |
| train | 0.50 | fitting all models, fitting all transforms | — |
| calibration | 0.25 | conformal quantiles only | any fitting or tuning |
| test | 0.25 | final evaluation, read once | anything that influences a choice |

**Why a separate calibration set.** Split conformal prediction needs nonconformity scores from
data the model has never seen. Scoring on training residuals gives optimistically small scores
and intervals that undercover. The calibration set is what buys the finite-sample guarantee,
and it is spent the moment it influences any modelling decision.

### Running stage 02

Stage 02 is lightweight and may run on the login node. While stage 01 is still a stub,
explicitly select the synthetic fixture:

```bash
make test  # generates the synthetic fixture if missing
make split SPLIT_ARGS="--synthetic --allow-dirty"
```

This writes `data/processed/splits.json` and
`outputs/manifests/02_make_splits_<UTC timestamp>.json`. These are **synthetic development
artefacts, never thesis results**. `--allow-dirty` explicitly permits the current uncommitted
setup; the manifest records `git_dirty: true` and `git_sha: null` until the first commit.
Without that flag, a dirty tree or an unborn HEAD fails before writing anything.

Once stage 01 exists, `make split` expects `<paths.processed>/cohort.csv`. An explicit
`--input /path/to/cohort.csv` is also supported. Inputs always go through `load_cohort` and
schema validation. There is no automatic fallback to synthetic data. Relative paths in the
configuration are relative to the repository; relative `--input` paths are relative to the
working directory.

The split file's version-1 contract is a JSON object with `schema_version`, `seed`,
`fractions` (train/calibration/test), `stratify_by`, `cohort_sha256`, `input_kind`, and `splits`.
`splits` contains exactly `train`, `calibration` and `test`, each a sorted list of patient-ID
strings. Both feature sets read these same lists. IDs that look numeric retain leading zeros.

Partition sizes use largest-remainder rounding of the configured fractions; ties follow
train, calibration, test order. Event counts are allocated proportionally to those sizes,
subject to including both events and censored rows in each partition. Infeasible tiny cohorts
fail explicitly. Each event stratum is sorted by ID and shuffled with
`np.random.default_rng(config.seed)`, so assignments do not depend on input row order.

Compatible reruns validate and reuse the existing file **without resampling or rewriting**.
They write a fresh manifest, while the split bytes and modification time stay unchanged.
Changed seed, fractions, stratification, source kind, schema version or cohort bytes cause an
error, as does a malformed or overlapping partition. Other configuration changes retain the
split and receive a new config hash in the run manifest. The cohort checksum covers the entire
input file, so even reordering CSV rows makes an existing file incompatible; newly generated
assignments remain invariant to row order.

There is no overwrite flag. Before transitioning from the synthetic check to real data,
archive the synthetic split and its manifests, or select separate processed/output directories
in a configuration file. Changing a real split invalidates every downstream fit and calibration
and needs an explicit experiment decision.

Manifests contain the Git SHA and dirty flag, SHA256 of the config/input/split files, resolved
package versions, seed, UTC timestamp, SLURM job ID (or null), elapsed wall time, source kind,
partition sizes and whether the frozen split was reused. They contain no patient IDs.

## The experiment matrix

Three models x two feature sets = six fitted models, all on the identical split.

| | clinical | clinical + PAM50 |
| --- | --- | --- |
| Cox PH | | |
| Bayesian Weibull AFT | | |
| Random Survival Forest | | |

RSF is deliberately left untuned — the brief rules out hyperparameter optimisation. State this
whenever RSF is compared against the other two.

## Research questions and the artefact that answers each

| RQ | Question | Artefact |
| --- | --- | --- |
| 1 | How well do clinical characteristics predict survival? | `outputs/metrics/performance.json`, clinical arm |
| 2 | Does PAM50 improve the prognosis? | `outputs/metrics/rq2_feature_set_comparison.json` |
| 3 | How do classical / Bayesian / ML models differ? | performance table across all six cells |
| 4 | What uncertainty does the Bayesian model give? | `outputs/metrics/bayes_diagnostics_*.json`, posterior predictive intervals |
| 5 | Do conformal intervals achieve 80/90/95%? | `outputs/metrics/conformal_coverage.json` |
| 6 | Is the best-performing model also the best-calibrated? | joint table: C-index vs coverage vs width |

RQ6 is the thesis' point: discrimination and trustworthy uncertainty are different properties
and need not coincide.

## TODO once results exist

- [ ] Record the realised split sizes and event rates per set.
- [ ] Record the conclusion for each RQ here, including null results.
