---
name: distill-paper-note
description: Turn a full-text paper in references/literature/ into a short structured note capturing method, assumptions and relevance to this thesis. Use when asked to summarise, review, read or take notes on a paper, to explain what a reference says, or to find which references support a claim.
---

# Distil a paper into a note

`references/literature/*.md` are **full-text PDF conversions**, 130-3,000 lines each, with OCR
artefacts (stray `<sup>` tags, `![](_page_3.jpeg)`, page numbers mid-sentence, and in three
files no headings at all). They are a corpus, not notes. This skill produces the note.

## 1. Navigate, do not read end to end

Loading a 3,000-line paper to answer one question wastes most of it. Instead:

1. Read the first ~80 lines for title, authors, abstract.
2. `grep -n` for the section you need — `"^## "` to list the paper's own headings, then jump.
3. Pull the method and assumptions; skip related-work and funding sections.
4. Expect OCR damage in `kumar_2019`, `unal_2026` and `yang_2024` — they have no headings and
   word-joining errors (`RESEARCHARTICLE`). Grep on partial words there.

## 2. Write the note

Use the template in `docs/literature-note-template.md`, saved to
`references/notes/<same_stem_as_the_paper>.md` so note and source line up.

Fill, in this order:
1. **Citation** — authors, year, venue.
2. **Problem** — one sentence.
3. **Method** — enough to reimplement or to say why you cannot.
4. **Assumptions** — the load-bearing part. For a conformal paper: what exchangeability is
   assumed, how censoring is handled, whether coverage is marginal or conditional. For a
   survival paper: censoring mechanism, proportional hazards or not, competing risks.
5. **Relevance here** — tie to a specific research question (RQ1-6) or say it is background.
6. **Whether it is usable** — is there code? Which language? `cfsurvival` and Gui's adaptive
   cut-offs are **R**, so using them means reimplementing, not importing. Record that.
7. **Quotes with line numbers**, so a claim can be traced back without rereading.

## 3. Be accurate about what the paper does and does not claim

- Distinguish what the authors proved from what they observed empirically.
- Record the conditions of their guarantee. "Valid coverage" almost always means marginal
  coverage under exchangeability — not conditional, and not under distribution shift.
- Note the evaluation setting. A method validated on uncensored simulation may not transfer
  to METABRIC's censoring.
- Never soften a limitation to make the paper look more applicable than it is.

## 4. Link it back

If the note changes a design decision, update the relevant `docs/` file and say which decision
moved. A note nothing points at will not be read again.
