#!/usr/bin/env bash
# PreToolUse(Bash): refuse heavy compute on the shared SLURM login node.
#
# openlab-slurm.hpc.local is a head node with 8 cores / 15 GB shared by the
# whole group. Four NUTS chains or a 500-tree forest degrade it for everyone.
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

payload="$(cat)"
cmd="$(printf '%s' "$payload" | extract_cmd)"
[ -n "$cmd" ] || exit 0

# Inside a SLURM allocation the cores are genuinely ours.
[ -n "${SLURM_JOB_ID:-}" ] && exit 0

case "$(hostname)" in
  *slurm*|*login*|*head*) ;;
  *) exit 0 ;;
esac

scan="$(printf '%s' "$cmd" | strip_heredocs)"

# Already a submission -- exactly the behaviour we want.
printf '%s' "$scan" | grep -Eq '(^|[;&|][[:space:]]*)[[:space:]]*(sbatch|srun|salloc)[[:space:]]' && exit 0

SUGGEST="Submit it instead: 'sbatch scripts/slurm/run_stage.sh <NN>', then poll with 'squeue -u \$USER'. See docs/slurm-guide.md."

if printf '%s' "$scan" | grep -Eq 'make[[:space:]]+(all|fit|conformal|eval|test-all)([[:space:]]|$)'; then
  deny "Blocked: that make target is compute, and this is the shared SLURM login node (8 cores / 15 GB). $SUGGEST"
fi

# Only fires when the stage is being EXECUTED -- reading or grepping one is fine.
if printf '%s' "$scan" | grep -Eq '(^|[;&|[:space:]])(python3?|ipython|\./scripts/)' \
   && printf '%s' "$scan" | grep -Eq 'scripts/0[345]_[a-z_]*\.py'; then
  deny "Blocked: pipeline stages 03-05 (fit / conformal / evaluate) are compute and must not run on the login node. $SUGGEST"
fi

# Only fires when Python is actually being invoked -- not when a file that
# mentions these names is merely being written or read.
if printf '%s' "$scan" | grep -Eq '(^|[;&|[:space:]])(python3?|ipython|jupyter)([[:space:]]|$)' \
   && printf '%s' "$scan" | grep -Eq '(pm\.sample|import[[:space:]]+pymc|nuts_sampler|RandomSurvivalForest\()'; then
  deny "Blocked: MCMC sampling and forest fitting must not run on the shared login node. $SUGGEST"
fi

if printf '%s' "$scan" | grep -Eq 'pytest[^|;]*(-m[[:space:]]*["'"'"']?slow([[:space:]"'"'"']|$)|--runslow)'; then
  deny "Blocked: the 'slow' tests sample real chains. Run them under sbatch, not on the login node. $SUGGEST"
fi

exit 0
