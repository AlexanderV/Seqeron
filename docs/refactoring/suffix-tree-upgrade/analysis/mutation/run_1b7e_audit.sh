#!/bin/bash
# Runs the hand-made mutants of the production code of 1b.7e in N separate git worktrees at the same time (the working checkout is never mutated); each worktree takes whole source files (the runner mutates one file at a
# time and restores it), the files are dealt to the workers by size. Every build and test run has a watchdog (run_1b7e_mutants.py kills the whole process group; a hang counts as killed). Results:
# results_1b7e_<round><file>.tsv next to this script. Do not run a timing measurement at the same time.
#   usage: run_1b7e_audit.sh <worktree prefix> <workers> [round prefix: "", "r2_", "r3_"]
set -u
PREFIX=${1:-/tmp/1b7e-audit}
N=${2:-3}
R=${3:-}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
for k in $(seq 1 "$N"); do
  W=$PREFIX-$k
  [ -d "$W" ] || git -C "$ROOT" worktree add --detach "$W" HEAD > /dev/null
  git -C "$W" checkout --detach -q "$(git -C "$ROOT" rev-parse HEAD)"
  mkdir -p "$W/tests/SuffixTree/SuffixTree.Esa.Slim"
  cp "$HERE/SuffixTree.Esa.Slim1b7e.csproj.template" "$W/tests/SuffixTree/SuffixTree.Esa.Slim/SuffixTree.Esa.Slim.csproj"
done
python3 -I - "$HERE" "$R" "$N" <<'PYEOF' > "$HERE/.deal_1b7e.tmp"
import glob, json, os, sys
here, rnd, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
files = []
for f in glob.glob(os.path.join(here, '1b7e_%s*_listed.json' % rnd)):
    base = os.path.basename(f)
    if not rnd and base.startswith('1b7e_r'):
        continue
    files.append((len(json.load(open(f))), f))
loads = [0] * n
deal = [[] for _ in range(n)]
for size, f in sorted(files, reverse=True):
    k = loads.index(min(loads))
    deal[k].append(f)
    loads[k] += size
for k in range(n):
    print(k + 1, ' '.join(deal[k]))
PYEOF
while read -r K FILES; do
  (
    for L in $FILES; do
      NAME=$(basename "$L" _listed.json)
      python3 -I "$HERE/run_1b7e_mutants.py" "$PREFIX-$K" "$L" "$HERE/results_$NAME.tsv"
    done
    echo "worker $K done" > "$HERE/.worker_1b7e_${R}$K.done"
  ) &
done < "$HERE/.deal_1b7e.tmp"
wait
rm -f "$HERE/.deal_1b7e.tmp"
echo "audit-done $R" > "$HERE/1b7e_${R}audit.done"
