# conformal-brca-prognosis

Trustworthy prognostic modelling for breast cancer using Bayesian survival analysis and
conformal prediction, on the METABRIC cohort.

Three models (Cox proportional hazards, Bayesian Weibull AFT, Random Survival Forest) are
compared across two feature sets (clinical, clinical + PAM50) on a fixed
train / calibration / test split, then wrapped in split conformal prediction to test whether
the declared 80% / 90% / 95% intervals achieve their stated coverage.

```bash
make setup    # uv sync
make test     # fast tests
make help     # all targets
```

Start with [CLAUDE.md](CLAUDE.md) for conventions and [docs/experimental-design.md](docs/experimental-design.md)
for the study design. The project brief is in [instructions.md](instructions.md).

> METABRIC is access-controlled (EGA `EGAS00000000083`). No patient data is in this repository,
> and none may be committed to it.
