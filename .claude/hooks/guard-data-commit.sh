#!/usr/bin/env bash
# PreToolUse(Bash): refuse to stage or commit patient data or generated outputs.
#
# METABRIC is access-controlled (EGA EGAS00000000083 / Synapse syn1688369) and
# git history is permanent: a patient row committed once is not undone by a
# later `git rm`. .gitignore is the first line of defence; this is the second.
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

payload="$(cat)"
cmd="$(printf '%s' "$payload" | extract_cmd)"
[ -n "$cmd" ] || exit 0

cmd="$(printf '%s' "$cmd" | strip_heredocs)"

# Only inspect commands that actually stage or commit.
printf '%s' "$cmd" | grep -Eq 'git[[:space:]]+(add|commit|stash[[:space:]]+push)' || exit 0

if printf '%s' "$cmd" | grep -Eq 'git[[:space:]]+add[[:space:]]+(-A|--all|-u|\.)([[:space:]]|$)'; then
  deny "Blanket 'git add' is blocked here. data/ holds access-controlled METABRIC records and git history is permanent. Stage by explicit path instead, e.g. 'git add src/brca/conformal/split.py'."
fi

if printf '%s' "$cmd" | grep -Eq 'git[[:space:]]+add[[:space:]]+.*(-f|--force)([[:space:]]|$)'; then
  deny "'git add --force' is blocked: it bypasses the .gitignore entries that keep data/ and outputs/ out of history. If a file genuinely belongs in git, move it out of those directories first."
fi

if printf '%s' "$cmd" | grep -Eq '(^|[[:space:]"'"'"'/])(data|outputs)/'; then
  deny "Blocked: this command touches data/ or outputs/. Neither may enter git -- data/ is access-controlled patient data, outputs/ is fully regenerable from code + config + seed. See .claude/rules/data-and-privacy.md."
fi

if printf '%s' "$cmd" | grep -Eq '\.(csv|tsv|parquet|nc|feather|xlsx|h5|pkl)([[:space:]"'"'"']|$)'; then
  deny "Blocked: staging a data file. If it is genuinely synthetic and not patient-derived, it belongs under tests/fixtures/ -- stage it from there by explicit path."
fi

# Notebooks: block only when cell outputs are actually present, so the
# suggested fix (strip, then stage) is reachable rather than a dead end.
if printf '%s' "$cmd" | grep -Eq '\.ipynb'; then
  for nb in $(printf '%s' "$cmd" | tr ' ' '\n' | grep -E '\.ipynb' | tr -d '"'"'"'"'); do
    if [ -f "$nb" ] && grep -q '"output_type"' "$nb" 2>/dev/null; then
      deny "Blocked: notebook '$nb' still contains cell outputs, which can embed patient-level rows. Run 'uv run nbstripout \"$nb\"' first, then stage it."
    fi
  done
fi

exit 0
