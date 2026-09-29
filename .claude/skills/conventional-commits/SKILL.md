---
name: conventional-commits
description: Choose the correct Conventional Commits type and draft a properly formatted commit message/subject line. Use when writing or suggesting a commit message, deciding between commit types like feat/fix/chore/refactor/style/build/ci, or when the user asks how a change should be committed. Does not run git commit itself — this repo's CLAUDE.md reserves committing/pushing for the user.
---

# Conventional Commits Cheatsheet

Use this to pick the right type and draft a commit message. Never run `git commit`
or `git push` — only draft the message text for the user to use themselves.

## Format

```
<type>[optional scope]: <description>
```

## Decision guide

Ask in this order:

1. **Reverting a prior commit?** → `revert` (reference the original commit's subject/hash)
2. **Only whitespace/formatting, no logic change?** → `style`
3. **New feature or capability the user didn't have before?** → `feat`
4. **Fixing broken/incorrect behavior?** → `fix`
5. **Only test files changed?** → `test`
6. **Only docs/comments/README/CHANGELOG?** → `docs`
7. **Restructuring code with no behavior change (extract function, rename, simplify, add null-safety)?** → `refactor`
8. **Measurable speed/memory improvement?** → `perf`
9. **CI/CD config (GitHub Actions, GitLab CI, Jenkins)?** → `ci`
10. **Build tooling or production dependency (Webpack, Docker, npm/Maven prod deps)?** → `build`
11. **Everything else (chores, dev-dependency bumps, .gitignore, scripts, dev-only suppression annotations)?** → `chore`

## Types reference

| Type | Use for | Examples |
|---|---|---|
| `feat` | New feature/functionality/behavior/config option | `feat(api): add support for pagination in user endpoint` |
| `fix` | Correcting broken/incorrect behavior in production code | `fix: null pointer handling` |
| `perf` | Measurable performance improvement | `perf: reduce number of redundant API calls` |
| `refactor` | Structural code change, no behavior change | `refactor: extract utility functions for data validation` |
| `style` | Cosmetic-only change (formatting, whitespace, indentation) | `style: reformat code with ESLint rules` |
| `test` | Adding/fixing/improving tests | `test(auth): improve token validation tests` |
| `docs` | Documentation, comments, API descriptions | `docs: update README to include installation steps` |
| `build` | Build process or production dependency changes | `build(deps): update express to v4.18.1` |
| `ci` | CI/CD configuration or workflow changes | `ci: add code quality checks in GitHub Actions` |
| `chore` | Admin/maintenance not affecting production code | `chore(deps): update eslint to v8.14.0` |
| `revert` | Rolling back a previous commit | `revert: "feat: add social login feature"` |

## Key distinctions (the grey areas)

- **`refactor` vs `style`**: `refactor` changes code structure/logic without changing behavior (extract function, add `Optional`/null-safety, simplify a loop). `style` changes only appearance (formatting, spacing) — zero structural impact.
- **`build` vs `chore`**: `build` is for changes that affect the production build or **production** dependencies. `chore` is for dev-dependency bumps, `.gitignore`, local scripts, and other things that don't touch what ships.
- **`fix` vs `refactor`**: `fix` corrects something that's actually broken/incorrect. `refactor` improves code that already works correctly.

## Workflow

1. Read the diff/change being committed (don't guess from the request alone).
2. Walk the decision guide above to pick the type.
3. Add a scope in parentheses when the change is localized to one area (e.g. `feat(ui): ...`), omit it when the change is broad or repo-wide.
4. Write the description in imperative mood, lowercase, no trailing period (e.g. `add`, not `added`/`adds`).
5. Present the drafted message to the user — do not run `git commit` yourself.

Source: adapted from [Nicola Bava's Conventional Commit Types cheatsheet](https://bavaga.dev) (Jan 2025), cross-referencing commit history from Angular, Electron, Jenkins X, and yargs.
