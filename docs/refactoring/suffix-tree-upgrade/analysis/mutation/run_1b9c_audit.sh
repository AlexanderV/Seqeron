#!/bin/bash
# Runs the lists 1b9c_<d|e>_<file>_listed.json (the mutants of 1b.7d and 1b.7e that the audits did not run, make_1b9c_lists.py) in N separate git worktrees at the same time (the working checkout is never mutated); each worktree takes
# whole lists (one source file each), dealt to the workers by size. The runner (run_1b9c_mutants.py) appends to results_1b9c_<d|e>_<file>.tsv after every mutant and skips the mutants that a result file already holds, so the
# script can be started again after a restart of the container and goes on where it stopped. Do not run a timing measurement at the same time.
#   usage: run_1b9c_audit.sh <worktree prefix> <workers> [list glob, default "1b9c_[de]_*_listed.json"]   (stage 2: AUDIT_COMMIT=<head of the audit> MUTANT_STAGES=12 MUTANT_SKIP_LISTED=1 MUTANT_TEST_LIMIT=1500 run_1b9c_audit.sh /tmp/1b9c-w 3 "1b9c_s2_*_listed.json")
# environment: MUTANT_STAGES (default "1"), MUTANT_BUILD_LIMIT, MUTANT_TEST_LIMIT, MUTANT_TMPDIR_BASE (the temporary files of the tests of worker K go to <base>K, default /dev/shm/1b9c-t: a memory file system, 3x faster) (see the runner)
set -u
PREFIX=${1:-/tmp/1b9c-w}
N=${2:-2}
GLOB=${3:-1b9c_[de]_*_listed.json}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
unset DOTNET_GCHeapHardLimit
HEADC=${AUDIT_COMMIT:-$(git -C "$ROOT" rev-parse HEAD)}  # AUDIT_COMMIT: the commit that the worktrees check out (stage 2 runs at the head of the audit, before the new tests)
for k in $(seq 1 "$N"); do
  W=$PREFIX$k
  [ -d "$W" ] || git -C "$ROOT" worktree add --detach "$W" HEAD > /dev/null
  git -C "$W" checkout --detach -q "$HEADC"
done
python3 -I - "$HERE" "$GLOB" "$N" <<'PYEOF' > "$HERE/.deal_1b9c.tmp"
import glob, json, os, sys
here, pat, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
files = []
for f in glob.glob(os.path.join(here, pat)):
    name = os.path.basename(f)[:-len('_listed.json')]
    res = os.path.join(here, 'results_%s.tsv' % name)
    done = set(l.split('\t')[0] for l in open(res)) if os.path.exists(res) else set()
    todo = [m for m in json.load(open(f)) if m['id'] not in done]
    if todo:
        files.append((len(todo), f))
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
      MUTANT_TMPDIR=${MUTANT_TMPDIR_BASE:-/dev/shm/1b9c-t}$K python3 -I "$HERE/run_1b9c_mutants.py" "$PREFIX$K" "$L" "$HERE/results_$NAME.tsv"
    done
    echo "worker $K done" > "$HERE/.worker_1b9c_$K.done"
  ) &
done < "$HERE/.deal_1b9c.tmp"
wait
rm -f "$HERE/.deal_1b9c.tmp"
echo "audit-done" > "$HERE/1b9c_audit.done"
