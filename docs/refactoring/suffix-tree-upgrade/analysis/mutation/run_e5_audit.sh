#!/bin/bash
# Runs the hand-made mutants of the new E5 code, one production file after the other, in a separate git worktree (E5_ROOT) so that the working checkout is never mutated.
# Every build and test run has a watchdog (run_e5_mutants.py kills the whole process group; a hang counts as killed). Results: results_e5_<file>.tsv next to this script.
#   usage: run_e5_audit.sh <worktree> [round tag]
set -u
W=${1:-/tmp/stbench/wt_e5}
TAG=${2:-}
HERE=$(cd "$(dirname "$0")" && pwd)
export E5_ROOT=$W PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
LIST="Ph5/SpillStack.cs:spillstack Ph5/AlignedWindowWriter.cs:windowwriter Ph5/CompactLcpFile.cs:compactlcpfile Ph5/ChildScan.cs:childscan Ph5/NsvScan.cs:nsvscan Ph5/ExternalChild.cs:externalchild Ph5/NextSymScan.cs:nextsym Experiments/E5.cs:e5 Harness/DirtySampler.cs:dirtysampler"
for item in $LIST; do
  F=${item%%:*}; N=${item##*:}
  python3 -I "$HERE/run_e5_mutants.py" "apps/SuffixTree.EmProto/$F" "$HERE/e5_${TAG}${N}_mutants.json" "$HERE/results_e5_${TAG}${N}.tsv"
done
echo audit-done > "$HERE/e5_${TAG}audit.done"
