#!/usr/bin/env bash
# Stop: refuse to declare work finished while fast tests fail or sensitive
# files are staged. This is the "must happen regardless of what Claude
# decides" gate -- it runs whether or not the tests were remembered.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO" || exit 0

payload="$(cat)"
# Do not re-block a continuation this hook itself triggered.
printf '%s' "$payload" | jq -e '.stop_hook_active == true' >/dev/null 2>&1 && exit 0

block() {
  jq -n --arg r "$1" '{decision:"block", reason:$r}'
  exit 0
}

# 1. Nothing sensitive staged.
staged="$(git diff --cached --name-only 2>/dev/null || true)"
if [ -n "$staged" ]; then
  bad="$(printf '%s\n' "$staged" | grep -E '^(data|outputs)/|\.(csv|tsv|parquet|nc|feather|h5|pkl)$' || true)"
  if [ -n "$bad" ]; then
    block "These staged files must never be committed: $(printf '%s' "$bad" | tr '\n' ' '). Unstage them with 'git restore --staged <path>' before finishing. See .claude/rules/data-and-privacy.md."
  fi
  for nb in $(printf '%s\n' "$staged" | grep -E '\.ipynb$' || true); do
    if [ -f "$nb" ] && grep -q '"output_type"' "$nb" 2>/dev/null; then
      block "Staged notebook '$nb' still has cell outputs, which can embed patient-level rows. Run 'uv run nbstripout \"$nb\"' and re-stage."
    fi
  done
fi

# 2. STATUS.md must reflect this session's work. It is the only memory the
#    other agent gets, so a change it does not record is effectively lost.
if [ -f STATUS.md ]; then
  stale="$(git status --porcelain --untracked-files=all 2>/dev/null \
    | sed -E 's/^.. //; s/^.* -> //' \
    | grep -vE '^(STATUS\.md|data/|outputs/|logs/|\.venv/)|(^|/)__pycache__/|\.pyc$' \
    | while IFS= read -r f; do [ -f "$f" ] && [ "$f" -nt STATUS.md ] && echo "$f"; done \
    | head -5 | tr '\n' ' ')"
  if [ -n "$stale" ]; then
    block "Files changed after STATUS.md was last updated: ${stale}. Update STATUS.md (milestones, decisions, a Handoff log entry, release your claim) per the status-handoff skill -- or, if nothing here changes project status, say so explicitly."
  fi
fi

# 3. Fast tests must pass. Skipped before the env exists, so `make setup` is
#    reachable on a fresh clone.
if [ -d .venv ] && [ -d tests ]; then
  if ! out="$(uv run --quiet pytest -m "not slow" -q 2>&1)"; then
    block "Fast tests are failing, so this work is not done. Fix them or say explicitly why they are expected to fail. Tail of output:
$(printf '%s' "$out" | tail -25)"
  fi
fi

exit 0
