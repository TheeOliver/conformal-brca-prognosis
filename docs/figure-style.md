# Figure and table style

> Pulled in by: `export-thesis-figure` skill; `.claude/rules/figures-and-tables.md`.

Target is **print**: matplotlib to vector PDF, dropped into a thesis document written
elsewhere. Screen conventions (hover, dark mode, tooltips) do not apply; print constraints
(greyscale photocopying, small reproduction size, fixed text width) do.

## Output

| Setting | Value |
| --- | --- |
| Format | vector **PDF** for the document, plus a PNG preview at the same basename |
| Width | one text width, set once here — never scale in the document |
| Font | match the thesis body font; label size >= 8pt at final printed size |
| Title | **none baked in** — the caption lives in the thesis where it can be edited |
| Name | by the claim: `rq5_coverage_vs_nominal.pdf`, not `plot3.pdf` |

## Categorical palette — the five PAM50 subtypes

These are the first five slots of a pre-validated colourblind-safe categorical order. Assign
in **fixed order**; never cycle, never recolour when a filter drops a subtype.

| Slot | Subtype | Hex | Linestyle |
| --- | --- | --- | --- |
| 1 | LumA | `#2a78d6` | solid |
| 2 | LumB | `#eb6834` | dashed |
| 3 | Her2 | `#1baf7a` | dash-dot |
| 4 | Basal | `#eda100` | dotted |
| 5 | Normal | `#e87ba4` | long-dash |

Two constraints came with these values and both apply here:

- **Colour is never the only channel.** Pair every hue with the linestyle above. This is
  required for print anyway — a photocopied figure is greyscale.
- **Relief rule.** Aqua, yellow and magenta sit below 3:1 contrast on a light surface, so
  those series need **visible direct labels** (or the accompanying table), not a legend alone.

Grey (`#767676`) is reserved for censoring marks, reference lines and at-risk annotations.
Never use it for a subtype.

> Note: the palette validator could not be run in this environment (no `node`), so these are
> the reference instance's published values rather than a locally re-validated set. If the
> palette is ever changed, run the `validate_palette.js` script bundled with the
> dataviz skill (external to this repo) before adopting it — do not eyeball colourblind safety.
>
> Field convention: some breast-cancer literature uses its own PAM50 colours. If you adopt
> this palette instead, say so once in the caption so an examiner is not misled.

## Marks

| Element | Spec |
| --- | --- |
| Lines | 2pt, round join and cap |
| Markers / end dots | >= 8pt diameter, filled in the series colour |
| Gridlines and axes | hairline, **solid** (never dashed), one step off the surface, recessive |
| Overlapping marks | 2pt ring in the surface colour so they separate |
| Labels | direct labels before a legend; a legend before a second axis |

Axis labels carry units — "Time (months)", not "time". Round tick values to clean numbers.

## Figure recipes

**Kaplan-Meier.** Step curves (`drawstyle="steps-post"`), censoring ticks in grey, an at-risk
table beneath aligned to the x-ticks, and a shaded confidence band at low alpha. Direct-label
each curve at its right-hand end rather than relying on the legend.

**Coverage vs nominal.** Nominal on x, empirical on y, with the **y = x reference line in
grey** so over- and under-coverage read at a glance. Binomial confidence intervals as error
bars — a point at 0.88 on 200 patients is not distinguishable from 0.90, and the figure should
show that rather than implying a miss.

**Coverage paired with width.** **Two stacked panels sharing the x-axis** — never a dual-axis
chart. Two y-scales on one plot is the most common chart error there is, and here it would
actively obscure the validity/efficiency trade-off the thesis is about. Coverage on top, width
(median with IQR band) below.

**Model comparison.** Dot-and-interval, sorted by the point estimate, with the bootstrap CI as
the interval. Not bars — bars imply a zero baseline that a C-index does not have.

## Tables

LaTeX `booktabs` into `outputs/tables/*.tex`. No vertical rules. Decimal-aligned. Uncertainty
beside every estimate. Consistent precision down a column — never 0.7123 beside 0.71.

## Before handing over

1. Legible at final print size, not just on screen.
2. Readable in greyscale.
3. Every number traces to a metrics JSON (`.claude/rules/figures-and-tables.md`).
4. Suggested caption text uses associative, non-causal language.
