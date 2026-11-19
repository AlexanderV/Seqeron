#!/usr/bin/env bash
# Audit of the harness of work package 1b.8 (outside the shipped code): every mutant of the lists 1b8_*_listed.json is built in a separate git worktree and run against the budgeted harness (FUZZ_BUDGET_S seconds, default 60).
# Two worktrees at the head of the checkout (commit the harness first) run side by side (FUZZ_THREADS = 2 workers of the harness each); the lists are dealt out alternately. Results: results_1b8_<list>.tsv (+ one .rep.log per mutant with the failure that killed it).
# usage: run_1b8_audit.sh [prefix of the lists, default 1b8_]       environment: AUDIT_WORK (worktree directory, default /tmp/claude-0/wt1b8), ROUND (a prefix like r2 for a re-run: lists 1b8_r2_*)
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE=$(pwd)
ROOT="$(cd ../../../../.. && pwd)"
WORK=${AUDIT_WORK:-/tmp/claude-0/wt1b8}
PREFIX=${1:-1b8_}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
mkdir -p "$WORK"
HEAD_COMMIT=$(git -C "$ROOT" rev-parse HEAD)
for w in a b; do
  if [ ! -d "$WORK/$w" ]; then git -C "$ROOT" worktree add --detach "$WORK/$w" "$HEAD_COMMIT" > /dev/null 2>&1 || { echo "worktree $w failed"; exit 1; }
  else git -C "$WORK/$w" checkout -q --detach "$HEAD_COMMIT" || { echo "checkout in $w failed"; exit 1; }; fi
done
lists=($(ls ${PREFIX}*_listed.json | sort))
i=0
declare -a A B
for f in "${lists[@]}"; do if [ $((i % 2)) -eq 0 ]; then A+=("$f"); else B+=("$f"); fi; i=$((i + 1)); done
run_set() { # worktree files...
  local wt=$1; shift
  for f in "$@"; do
    name=${f%_listed.json}
    python3 "$HERE/run_1b8_mutants.py" "$wt" "$HERE/$f" "$HERE/results_${name}.tsv" > "$HERE/results_${name}.log" 2>&1
  done
}
run_set "$WORK/a" "${A[@]}" &
run_set "$WORK/b" "${B[@]}" &
wait
touch "$HERE/${PREFIX}audit.done"
echo done
