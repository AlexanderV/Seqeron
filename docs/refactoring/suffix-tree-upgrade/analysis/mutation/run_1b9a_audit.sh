#!/usr/bin/env bash
# Audit of the code of work package 1b.9a (outside the shipped code): every mutant of the lists 1b9a_*_listed.json is built in a separate git worktree and run against the unit tests of the new code and, for a survivor,
# the budgeted harness of 1b.8 (FUZZ_BUDGET_S seconds, default 60). Two worktrees at the head of the checkout (commit first) run side by side; the lists are dealt out by size. Results: results_1b9a_<list>.tsv.
# usage: run_1b9a_audit.sh [prefix of the lists, default 1b9a_]       environment: AUDIT_WORK (worktree directory, default /tmp/wt1b9a), ROUND (a prefix like r2_ for a re-run: lists 1b9a_r2_*)
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE=$(pwd)
ROOT="$(cd ../../../../.. && pwd)"
WORK=${AUDIT_WORK:-/tmp/wt1b9a}
PREFIX=${1:-1b9a_}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
mkdir -p "$WORK"
HEAD_COMMIT=$(git -C "$ROOT" rev-parse HEAD)
for w in a b; do
  if [ ! -d "$WORK/$w" ]; then git -C "$ROOT" worktree add --detach "$WORK/$w" "$HEAD_COMMIT" > /dev/null 2>&1 || { echo "worktree $w failed"; exit 1; }
  else git -C "$WORK/$w" checkout -q --detach "$HEAD_COMMIT" || { echo "checkout in $w failed"; exit 1; }; fi
done
mapfile -t lists < <(ls ${PREFIX}*_listed.json | sort)
python3 -I - "$HERE" "${lists[@]}" <<'PYEOF' > "$HERE/.deal_1b9a.tmp"
import json, os, sys
here, files = sys.argv[1], sys.argv[2:]
sized = sorted(((len(json.load(open(os.path.join(here, f)))), f) for f in files), reverse=True)
loads = [0, 0]
deal = [[], []]
for size, f in sized:
    k = loads.index(min(loads))
    deal[k].append(f)
    loads[k] += size
for k in range(2):
    print('ab'[k], ' '.join(deal[k]))
PYEOF
run_set() { # worktree files...
  local wt=$1; shift
  for f in "$@"; do
    name=${f%_listed.json}
    python3 "$HERE/run_1b9a_mutants.py" "$wt" "$HERE/$f" "$HERE/results_${name}.tsv" > "$HERE/results_${name}.log" 2>&1
  done
}
while read -r W FILES; do
  run_set "$WORK/$W" $FILES &
done < "$HERE/.deal_1b9a.tmp"
wait
rm -f "$HERE/.deal_1b9a.tmp"
touch "$HERE/${PREFIX}audit.done"
echo done
