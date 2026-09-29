---
description: Split conformal prediction under censoring - finite-sample correction, calibration hygiene, honest reporting.
paths:
  - "src/brca/conformal/**"
  - "tests/test_conformal*.py"
---

# Conformal prediction

**This code is hand-written.** No Python library does conformal survival prediction under
censoring — the Candes and Gui reference implementations are R. Treat every line as unverified.

- **The quantile is `ceil((n + 1) * (1 - alpha)) / n`** over the `n` calibration scores — not
  `np.quantile(scores, 1 - alpha)`. The naive version undercovers by roughly `1/n`: invisible
  in a plot, wrong in the guarantee. When the index exceeds `n` the correct interval is
  **infinite** — return it and report it, never clip.
- `alpha` is miscoverage. A function called with `alpha=0.05` returns a **95%** interval.
- **Split (inductive) CP only.** Scores come from the calibration set under a model fit on
  train alone. Never fit on calibration, never tune on it, and never pick an alpha by looking
  at the resulting coverage.
- Censoring breaks plain exchangeability. State in the docstring which adaptation is used
  (IPCW-weighted scores, or a Candes-style lower predictive bound) and what it assumes.
- Every nonconformity score needs a unit test asserting that empirical coverage on synthetic
  exchangeable data lands within Monte-Carlo error of nominal. Coverage is the contract.
- Always return interval **width** with the interval: validity without efficiency is vacuous.
- Report the calibration set size `n` with every interval: the guarantee is finite-sample and
  meaningless without it.

Methods, assumptions and sources: `docs/conformal-methods.md`.
