#!/bin/bash
# Runs the hand-made mutants of the new E4 code, one production file after the other, in a separate git worktree (E4_ROOT) so that the working checkout is never mutated.
# Every build and test run has a watchdog (run_e4_mutants.py kills the whole process group; a hang counts as killed). Results: results_e4_<file>.tsv next to this script.
#   usage: run_e4_audit.sh <worktree> [first file index] [last file index]
set -u
W=${1:-/tmp/e4-audit}
HERE=$(cd "$(dirname "$0")" && pwd)
export E4_ROOT=$W PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
LIST="Ph4/LcpGeometry.cs:geometry Ph4/Spool.cs:spool Ph4/PlcpScan.cs:plcpscan Ph4/LcpOutput.cs:output Ph4/LcpBuild.cs:build Experiments/E4.cs:e4"
for item in $LIST; do
  F=${item%%:*}; N=${item##*:}
  python3 -I "$HERE/run_e4_mutants.py" "apps/SuffixTree.EmProto/$F" "$HERE/e4_${N}_mutants.json" "$HERE/results_e4_${N}.tsv"
done
echo audit-done > "$HERE/e4_audit.done"
