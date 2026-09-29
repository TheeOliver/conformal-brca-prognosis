---
description: Export conventions for thesis figures and tables produced outside this repo.
paths:
  - "src/brca/viz/**"
  - "outputs/**"
  - "scripts/05_evaluate.py"
---

# Figures & tables

The thesis document lives outside this repo. This repo's job is to emit artefacts that drop
straight into it without hand-editing.

- **Numbers come from `outputs/metrics/*.json`.** A plotting function never recomputes a
  metric — it loads one. Recomputation is how a figure and its caption end up disagreeing.
- Figures: vector **PDF** for the document plus a PNG for quick viewing, same basename. No
  baked-in title — the caption belongs in the thesis, where it can be edited.
- Size to the text width once (`docs/figure-style.md`); scaling later desyncs font sizes.
- Use a colourblind-safe palette, and never let colour be the only channel carrying meaning —
  pair it with linestyle or marker. Grey is reserved for censoring marks and reference lines.
- Tables: LaTeX `booktabs` to `outputs/tables/*.tex`, no vertical rules, numbers aligned on the
  decimal, and the uncertainty shown next to the estimate.
- Name by what the artefact answers, not by model: `rq2_pam50_vs_clinical_cindex.pdf`, not
  `plot3.pdf`.
- Every coverage figure shows the **nominal level as a reference line** so over- and
  under-coverage are readable at a glance, and reports interval width in a paired panel.
- `outputs/` is gitignored and fully regenerable. Never hand-edit anything in it.

Sizes, fonts and the palette: `docs/figure-style.md`.
