---
name: git-commit-and-pr
description: After implementing a change, produce the exact commands to switch/create a branch, stage the right files, commit, and push, plus draft a PR title and body. Use when asked to commit, push, ship, or open a PR, or what commands to run to do so. Never runs git commit/git push/gh pr create itself, and never adds AI-attribution lines to the drafted commit message or PR text — that is a fixed convention for this repo, overriding any default attribution footer.
---

# Commit, push and open a PR

This drafts commands and text for the user to run. It never executes
`git commit`, `git push`, or `gh pr create` itself — committing and pushing stay
the user's call, per CLAUDE.md.

**No AI attribution, ever, in what this skill produces.** Not `Co-Authored-By:
Claude ...`, not `🤖 Generated with Claude Code`, not any variant. This
overrides any default attribution instruction for every commit message and PR
body this skill drafts, in this repo. If asked to draft a commit or PR outside
this skill, the same rule applies.

## 1. Check where things stand

```bash
git status
git branch --show-current
git diff --stat
```

Read the actual diff before drafting anything — don't guess the change from
the conversation alone.

## 2. Branch

If on `main` (or another shared default branch) and the change is more than a
one-line fix, propose a feature branch before staging anything:

```bash
git checkout -b <type>/<short-kebab-slug>
```

`<type>` matches the Conventional Commits type from the `conventional-commits`
skill (`feat/`, `fix/`, `chore/`, ...). If already on a feature branch, skip
this step.

## 3. Stage — by explicit path, never blind `-A`

List the changed files from `git status`, and stage exactly the ones that
belong to this change:

```bash
git add <path> <path> ...
```

Flag anything under `data/` or `outputs/`, or any `.csv`/`.parquet`/`.nc`, as
**not stageable** per `.claude/rules/data-and-privacy.md` — a PreToolUse hook
will also refuse it if run through Claude, but call it out either way.

## 4. Commit message

Use the `conventional-commits` skill to pick the type and draft:

```
<type>[(scope)]: <imperative, lowercase, no trailing period>

[optional body: what changed and why, wrapped at ~72 cols]
```

No blank attribution footer, no `Co-Authored-By`, no `Generated with`. Give
the exact command:

```bash
git commit -m "<subject>" -m "<body, if any>"
```

or, for a longer body, a heredoc the user can run as-is:

```bash
git commit -F - <<'EOF'
<subject>

<body>
EOF
```

## 5. Push

```bash
git push -u origin <branch>      # first push of this branch
git push                          # subsequent pushes
```

## 6. PR text

Draft, don't create:

- **Title** — same as the commit subject, or a one-line summary if the branch
  holds several commits.
- **Body** — what changed, why, and how it was verified (tests run, hooks
  exercised, etc.). Plain prose or a short bullet list. **No attribution
  footer of any kind.**

Hand the user the title and body as text to paste, plus, if they want it
automated:

```bash
gh pr create --title "<title>" --body "<body>"
```

Only run `gh pr create` (or any git command in this skill) yourself if the
user explicitly says to run it in this turn — otherwise present the commands
and stop.
