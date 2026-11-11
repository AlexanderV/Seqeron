#!/bin/bash
# Runs the hand-made mutants of the new E3b code, one production file after the other, in a separate git worktree (E3B_ROOT) so that the working checkout is never mutated.
# Every build and test run has a watchdog (run_e3b_mutants.py kills the whole process group; a hang counts as killed). Results: results_e3b_<file>.tsv next to this script.
#   usage: run_e3b_audit.sh <worktree> [first file index] [last file index]
set -u
W=${1:-/tmp/e3b-audit}
HERE=$(cd "$(dirname "$0")" && pwd)
export E3B_ROOT=$W PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
LIST="Dcs/SampleGeometry.cs:samplegeometry Dcs/MappedWindowWriter.cs:mappedwindowwriter Dcs/Runs.cs:runs Dcs/SampleNamer.cs:samplenamer Dcs/NameSinks.cs:namesinks Dcs/SaInverter.cs:sainverter Dcs/ExternalDcs.cs:externaldcs Experiments/E3b.cs:e3b"
for item in $LIST; do
  F=${item%%:*}; N=${item##*:}
  python3 -I "$HERE/run_e3b_mutants.py" "apps/SuffixTree.EmProto/$F" "$HERE/e3b_${N}_mutants.json" "$HERE/results_e3b_${N}.tsv"
done
# the harness file is also used by the E0 tests: stage 2 is the whole EmProto test project
python3 -I "$HERE/run_e3b_mutants.py" apps/SuffixTree.EmProto/Harness/ArenaBudget.cs "$HERE/e3b_arenabudget_mutants.json" "$HERE/results_e3b_arenabudget.tsv" "" "FullyQualifiedName~EmProto.Tests"
echo audit-done > "$HERE/e3b_audit.done"
