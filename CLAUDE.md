# conformal-brca-prognosis — Claude Code

Project conventions are shared with Codex and live in `AGENTS.md`, imported here in full.
Edit conventions there, never here, so the two agents cannot drift apart.

@AGENTS.md

## Claude-specific

- A SessionStart hook injects the `STATUS.md` snapshot. Still read the full file before
  claiming work — the snapshot omits older handoff entries.
- Rules in `.claude/rules/` auto-load by path; the routing table in `AGENTS.md` exists for Codex.
- Skills are invocable directly. Reviewers in `.claude/agents/` run as subagents.
- Hooks (`.claude/settings.json`): auto-format on edit; deny data staging and login-node
  compute; the Stop gate blocks until fast tests pass and `STATUS.md` postdates your changes.
