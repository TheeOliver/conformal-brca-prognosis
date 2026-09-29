---
name: export-thesis-figure
description: Produce a publication-ready figure or LaTeX table from computed metrics, sized and styled for the thesis document. Use when asked to make, update or export a figure, plot, chart or table for the thesis, a chapter, a supervisor meeting, or a presentation.
---

# Export a thesis figure or table

The thesis document lives outside this repo. The job here is an artefact that drops in without
hand-editing. Read `docs/figure-style.md`; `.claude/rules/figures-and-tables.md` applies.

## 1. Locate the numbers — never recompute them

Find the metric in `outputs/metrics/*.json`. If it is not there, the fix is to run the stage
that produces it, **not** to compute it inline in the plotting code. A figure that recomputes
its own numbers will eventually disagree with the table beside it.

If the metrics file is stale relative to the code (check its manifest's git SHA), say so and
rerun rather than plotting old numbers.

## 2. Choose the form for the claim

| Claim | Form |
| --- | --- |
| model A ranks patients better than B | C-index with bootstrap CI, dot-and-interval |
| predictions are well calibrated | Brier over time, or calibration curve at a horizon |
| intervals achieve nominal coverage | coverage vs nominal, with the diagonal/reference line |
| intervals are informative | width distribution, paired with the coverage panel |
| survival differs by subtype | Kaplan-Meier with at-risk table and censoring ticks |

Coverage is **never** plotted alone — always pair it with width in the same figure, because a
valid-but-useless interval must be visible as such.

## 3. Build it

- Vector PDF plus a PNG preview, same basename, into `outputs/figures/`.
- No baked-in title — the caption belongs in the thesis where it can be edited.
- Size to text width once from `docs/figure-style.md`; never scale in the document.
- Colourblind-safe palette, and colour never the sole carrier of meaning — pair with
  linestyle or marker. Grey is reserved for censoring marks and reference lines.
- Axis labels carry units ("months", not "time").
- Name for the claim: `rq5_coverage_vs_nominal.pdf`, not `plot3.pdf`.

## 4. Tables

LaTeX `booktabs` into `outputs/tables/*.tex`. No vertical rules. Decimal-aligned. Uncertainty
next to every estimate. Consistent precision across a column — do not print 0.7123 beside 0.71.

## 5. Check before handing over

1. Legible at final print size, not just on screen.
2. Every number traces back to a metrics JSON.
3. Caption text you suggest uses associative, non-causal language.
4. State the filename and what it shows, so it can be dropped into the right chapter.
