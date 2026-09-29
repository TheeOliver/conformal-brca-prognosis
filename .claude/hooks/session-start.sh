#!/usr/bin/env bash
# SessionStart: put the shared STATUS.md snapshot in front of the agent before
# it does anything. Two agents alternate in this repo and neither remembers the
# other's session -- this is how each one learns where things stand.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
S="$REPO/STATUS.md"
cat >/dev/null   # drain the payload; nothing in it is needed

if [ ! -f "$S" ]; then
  echo "STATUS.md is missing. Recreate it from the status-handoff skill before starting work."
  exit 0
fi

echo "=== STATUS.md snapshot (read the full file before claiming work) ==="
# Everything above the handoff log, plus the two most recent handoff entries.
awk '
  /^## Handoff log/ { print; inlog = 1; next }
  inlog && /^### / { n++ }
  !inlog || n <= 2 { print }
' "$S"
echo "=== end snapshot ==="
exit 0
