#!/bin/bash
# Stage 2 of the audit of 1b.9c: the mutants that survived stage 1 run against the WHOLE default set of SuffixTree.Esa.Tests (every test that is not [Explicit]) in one worktree; the result files are
# results_1b9c_s2_<d|e>_<file>.tsv (appended, restartable). The lists are made again from the results of stage 1 first (make_1b9c_s2.py).
#   usage: run_1b9c_s2.sh <worktree> [slot]   (the slot names the directory of the temporary files)
set -u
W=${1:?worktree}
SLOT=${2:-s2}
HERE=$(cd "$(dirname "$0")" && pwd)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
unset DOTNET_GCHeapHardLimit
python3 -I "$HERE/make_1b9c_s2.py"
for L in "$HERE"/1b9c_s2_*_listed.json; do
  NAME=$(basename "$L" _listed.json)
  MUTANT_STAGES=2 MUTANT_TMPDIR=/dev/shm/1b9c-t$SLOT python3 -I "$HERE/run_1b9c_mutants.py" "$W" "$L" "$HERE/results_$NAME.tsv"
done
echo "s2-done" > "$HERE/1b9c_s2_$SLOT.done"
