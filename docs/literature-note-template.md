# Literature note template

> Pulled in by: `distill-paper-note` skill.

`references/literature/*.md` are full-text paper conversions, not notes. Notes are a **new
artefact**: write them to `references/notes/<same_stem_as_the_paper>.md` so note and source
line up by filename.

Naming follows the corpus convention, which is exceptionless across all 19 papers:
`{first_author_surname}_{year}_{snake_case_slug}.md`.

---

```markdown
---
key: candes_2023_conformalized_survival_analysis
authors: Candes, Lei, Ren
year: 2023
venue: JRSS-B
rq: [5]          # which research questions this bears on; [] for background
---

# Conformalized survival analysis

## Problem
One sentence: what could not be done before this paper.

## Method
Enough to reimplement, or an explicit note on why it cannot be reimplemented here.

## Assumptions
The load-bearing section. Be specific, and separate what is proved from what is observed.
- Conformal papers: what exchangeability is assumed; how censoring is handled; whether
  coverage is marginal or conditional; what happens under distribution shift.
- Survival papers: censoring mechanism; proportional hazards or not; competing risks;
  how the time grid was chosen.

## Relevance here
Tie to a specific RQ, or say plainly that it is background. A note with no tie will not be
read again.

## Usable?
Is there code? **Which language?** `cfsurvival` and Gui's adaptive cut-offs are R, so using
them means reimplementing, not importing. Record the repo and the language.

## Key quotes
> "..." (L1284)

Line numbers into the source file, so a claim is traceable without rereading.

## Limitations for this project
Where the method's conditions do not hold for METABRIC. Never soften these to make the paper
look more applicable than it is.
```

---

## Navigating the sources

The papers run 130-3,000 lines. Do not read one end to end to answer one question:

- first ~80 lines for title, authors, abstract;
- `grep -n "^## "` to list the paper's own headings, then jump;
- `kumar_2019`, `unal_2026` and `yang_2024` have **no headings** and OCR word-joining damage
  (`RESEARCHARTICLE`) — grep on partial words there.
