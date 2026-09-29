---
name: status-handoff
description: Read, claim and update STATUS.md, the shared progress ledger between Claude Code and Codex. Use at the start of every session, before starting any non-trivial task, and at the end of any session that changed files - or whenever asked "where are we", "what's next", "what did the other agent do", or to hand off work.
---

# STATUS.md handoff

Claude Code and Codex alternate in this repo and neither remembers the other's sessions.
`STATUS.md` is the only shared memory. A change it does not record is, for the other agent,
a change that never happened — and will likely be redone or undone.

## At session start

1. Read `STATUS.md` in full — not only the injected snapshot.
2. Read the newest *Handoff log* entry closely: its **Next** line is usually your task.
3. Check *Active claims*. If the task you were asked to do is claimed by the other agent:
   - claim dated today or yesterday → stop and tell the user; do not work around it;
   - older than ~3 days with no handoff since → ask the user whether it is abandoned.
4. `git branch --show-current` and `git status` — confirm the tree matches what the last
   handoff describes. If it does not (uncommitted changes nobody logged), tell the user before
   touching anything.

## Before non-trivial work — claim it

Add a row to *Active claims*:

```
| <task, one line> | Claude | <branch> | 2026-09-29 |
```

Use your real agent name (`Claude` or `Codex`) and today's date. Trivial edits (typo, one-line
fix) need no claim, but still need a handoff entry if they change files.

## At session end — record it

Update, in this order:

1. **Milestones** — tick what is actually done. Done means `make check` passes and the
   relevant test exists, not "code written".
2. **Decisions** — any choice the other agent must respect: an API, a prior, a threshold, a
   file layout, a rejected approach and why. Date and agent on each row.
3. **Blockers & open questions** — add new ones; strike through resolved ones with the
   resolution, do not delete them.
4. **Handoff log** — a new entry at the **top** of the log:

   ```
   ### 2026-09-29 — Claude — <short title>
   - **Changed:** files/areas touched, in one or two lines
   - **Verified:** commands actually run and their result (`make check` → 19 passed)
   - **Not done / caveats:** what is incomplete, stubbed, or unverified
   - **Next:** the single most useful next step, specific enough to start on cold
   ```

5. **Active claims** — remove your row.
6. **Snapshot** at the top — update *Phase*, *Last updated*, and *Next up*.

## Rules

- **Verified means run.** Never write "tests pass" without having run them this session. If
  you could not run something (e.g. it needs `sbatch`), write that under *Not done*.
- Never delete or rewrite the other agent's entries. If one is wrong, add a correction entry
  that says what was wrong.
- Keep entries short. `STATUS.md` is read at the start of every session by both agents; a
  bloated log costs both of them context. When the log passes ~15 entries, move the oldest
  into `docs/status-archive.md` in the same commit.
- No patient-level data in `STATUS.md` — aggregates only (convention 2).
- Numbers from the synthetic fixture are labelled as synthetic, always.
