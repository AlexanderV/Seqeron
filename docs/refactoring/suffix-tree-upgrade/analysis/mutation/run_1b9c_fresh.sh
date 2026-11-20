#!/bin/bash
# The fresh mutants of 1b.9c (make_1b9c_fresh.py: single-site mutants written after the new tests, over the code that the new tests touch) in one git worktree that holds the final code and the final tests:
# pass 1 runs the new test classes only (the question is whether the new tests kill them); every survivor of pass 1 runs the whole default set in pass 2 (an older test may kill it; a survivor of both is a gap).
# Results: results_1b9c_f_<file>.tsv (pass 1), results_1b9c_f2_<file>.tsv (pass 2), both appended and restartable.
#   usage: run_1b9c_fresh.sh <worktree>
set -u
W=${1:?worktree}
HERE=$(cd "$(dirname "$0")" && pwd)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
unset DOTNET_GCHeapHardLimit
for L in "$HERE"/1b9c_f_*_listed.json; do
  NAME=$(basename "$L" _listed.json)
  MUTANT_FILTER_ONLY='FullyQualifiedName~Gap1b9cTests' MUTANT_TMPDIR=/dev/shm/1b9c-tf python3 -I "$HERE/run_1b9c_mutants.py" "$W" "$L" "$HERE/results_$NAME.tsv"
done
python3 -I - "$HERE" <<'PYEOF'
import glob, json, os, sys
here = sys.argv[1]
for f in sorted(glob.glob(os.path.join(here, '1b9c_f_*_listed.json'))):
    name = os.path.basename(f)[:-len('_listed.json')]
    res = os.path.join(here, 'results_%s.tsv' % name)
    surv = set(l.split('\t')[0] for l in open(res) if l.split('\t')[2] == 'survived') if os.path.exists(res) else set()
    todo = [m for m in json.load(open(f)) if m['id'] in surv]
    if todo:
        json.dump(todo, open(os.path.join(here, '1b9c_f2_%s_listed.json' % name[len('1b9c_f_'):]), 'w'), indent=1)
PYEOF
for L in "$HERE"/1b9c_f2_*_listed.json; do
  [ -e "$L" ] || continue
  NAME=$(basename "$L" _listed.json)
  MUTANT_STAGES=2 MUTANT_TEST_LIMIT=1500 MUTANT_TMPDIR=/dev/shm/1b9c-tf python3 -I "$HERE/run_1b9c_mutants.py" "$W" "$L" "$HERE/results_$NAME.tsv"
done
echo "fresh-done" > "$HERE/1b9c_fresh.done"
