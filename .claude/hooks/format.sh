#!/usr/bin/env bash
# PostToolUse: keep edited Python formatted and notebooks output-free.
# Deliberately never blocks -- a formatter failure must not stop work.
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
payload="$(cat)"
cd "$REPO" || exit 0
[ -d .venv ] || exit 0   # nothing to run before `make setup`

# Claude edits one file per call; a Codex apply_patch can touch several.
printf '%s' "$payload" | extract_edited_files | while IFS= read -r file; do
  case "$file" in /*) ;; *) file="$REPO/$file" ;; esac
  [ -f "$file" ] || continue
  case "$file" in "$REPO"/*) ;; *) continue ;; esac   # only files in this project
  case "$file" in
    *.py)
      uv run --quiet ruff format "$file"      >/dev/null 2>&1
      uv run --quiet ruff check --fix "$file" >/dev/null 2>&1
      ;;
    *.ipynb)
      # Cell outputs can carry patient-level rows straight into git history.
      uv run --quiet nbstripout "$file"       >/dev/null 2>&1
      ;;
  esac
done
exit 0
