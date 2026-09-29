#!/usr/bin/env bash
# Shared helpers for the PreToolUse guards.

# Heredoc bodies are DATA, not commands. Writing a doc or a source file that
# merely mentions a guarded token must not be blocked -- only actually running
# the thing is. Strip heredoc bodies before pattern matching.
strip_heredocs() {
  python3 -c '
import re, sys
out, delim = [], None
for ln in sys.stdin.read().split("\n"):
    if delim is not None:
        if ln.strip() == delim:
            delim = None
        continue
    out.append(ln)
    m = re.search(r"<<-?\s*([\"\x27]?)([A-Za-z_][A-Za-z0-9_]*)\1", ln)
    if m:
        delim = m.group(2)
sys.stdout.write("\n".join(out))
'
}

# Emit a PreToolUse deny decision and stop.
deny() {
  jq -n --arg r "$1" '{hookSpecificOutput:{hookEventName:"PreToolUse",permissionDecision:"deny",permissionDecisionReason:$r}}'
  exit 0
}

# ---------------------------------------------------------------------------
# Payload extraction, shared by Claude Code and Codex.
#
# Claude:  Bash        -> tool_input.command (string)
#          Edit/Write  -> tool_input.file_path / notebook_path
# Codex:   exec_command-> tool_input.cmd (string)
#          exec        -> raw text
#          apply_patch -> raw patch text ("*** Update File: <path>")
# Unknown shapes yield nothing, so a guard fails open on tools it cannot read
# rather than blocking work it does not understand. The agent-independent
# layers (.githooks/pre-commit, src/brca/compute_guard.py) are the backstop.
# ---------------------------------------------------------------------------

# The shell command a tool call will run, or nothing for file-edit tools.
extract_cmd() {
  jq -r '
    (.tool_name // "") as $n
    | if ($n | test("^(apply_patch|Edit|Write|MultiEdit|NotebookEdit)$")) then ""
      else
        .tool_input as $t
        | ( if ($t|type) == "string" then $t
            elif ($t|type) == "object" then ($t.command // $t.cmd // $t.input // "")
            else "" end )
        | if type == "array" then join(" ") elif type == "string" then . else "" end
        | if startswith("*** Begin Patch") then "" else . end
      end'
}

# Files a tool call edited, one per line, absolute where resolvable.
extract_edited_files() {
  jq -r '
    (.cwd // "") as $cwd
    | [ (.tool_input | objects | (.file_path // .notebook_path // empty)),
        ( (.tool_input | if type == "string" then .
                         elif type == "object" then (.input // .patch // "")
                         else "" end)
          | if type == "string" then . else "" end
          | split("\n")[]
          | capture("^\\*\\*\\* (Add|Update) File: (?<p>.+)$")? | .p )
      ]
    | .[]
    | if startswith("/") or $cwd == "" then . else "\($cwd)/\(.)" end'
}
