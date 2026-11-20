#!/bin/bash
# Runs the lists 1b9c_c_<d|e>_<file>_listed.json (make_1b9c_closed.py: the mutants that a new test of 1b.9c closes) in one git worktree that holds the final code and the final tests: the only tests run are the new classes of this package
# (the kill must come from them: "the kill must reproduce"). Results: results_1b9c_c_<d|e>_<file>.tsv (appended, restartable).
#   usage: run_1b9c_closed.sh <worktree>
set -u
W=${1:?worktree}
HERE=$(cd "$(dirname "$0")" && pwd)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
unset DOTNET_GCHeapHardLimit
for L in "$HERE"/1b9c_c_*_listed.json; do
  NAME=$(basename "$L" _listed.json)
  MUTANT_FILTER_ONLY='FullyQualifiedName~SaGap1b9cTests|FullyQualifiedName~SesaGap1b9cTests' MUTANT_TMPDIR=/dev/shm/1b9c-tc python3 -I "$HERE/run_1b9c_mutants.py" "$W" "$L" "$HERE/results_$NAME.tsv"
done
echo "closed-done" > "$HERE/1b9c_closed.done"
