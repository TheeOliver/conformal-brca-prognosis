---
description: Language discipline - prognostic prediction is not causal inference, and coverage is not validity.
paths:
  - "docs/**"
  - "outputs/**"
  - "src/brca/evaluation/**"
  - "src/brca/viz/**"
---

# Scientific claims

This project is **prognostic**. Nothing in the design identifies a causal effect, and a strong
predictor need not affect anything. In comments, docstrings, docs, labels and captions:

| Do not write | Write instead |
| --- | --- |
| "the effect of tumour size on survival" | "tumour size is associated with survival" |
| "PAM50 improves patient outcomes" | "PAM50 improves predictive performance" |
| "grade 3 causes / drives / reduces risk" | "grade 3 is predictive of shorter survival" |
| "impact of nodal status" | "contribution of nodal status to the prediction" |

A Cox hazard ratio here is an association conditional on the other covariates, not an effect.

**Coverage is not validity by assertion.** Never claim an interval "is valid" or "guarantees
90%" from theory alone — the guarantee is conditional on exchangeability, which censoring
strains. State **measured** empirical coverage on test, with the calibration set size `n`.

**Report the negative result.** If clinical+PAM50 does not beat clinical, say so plainly and
prominently. A well-evidenced null answer to RQ2 is a real finding.
