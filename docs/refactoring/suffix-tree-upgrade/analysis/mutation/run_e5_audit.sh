#!/bin/bash
# Runs the hand-made mutants of the new E5 code, one production file after the other, in a separate git worktree (E5_ROOT) so that the working checkout is never mutated.
# Every build and test run has a watchdog (run_e5_mutants.py kills the whole process group; a hang counts as killed). Results: results_e5_<tag><file>.tsv next to this script.
#   usage: run_e5_audit.sh <worktree> <tag> <file:name[:slice]> ...      e.g. run_e5_audit.sh /tmp/wt "" Ph5/SpillStack.cs:spillstack Ph5/ChildScan.cs:childscan
#   two audits may run at the same time in two worktrees (4 vCPU): use disjoint files and different tags
set -u
W=$1; TAG=$2; shift 2
HERE=$(cd "$(dirname "$0")" && pwd)
export E5_ROOT=$W PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
for item in "$@"; do
  IFS=: read -r F N SL <<< "$item"
  python3 -I "$HERE/run_e5_mutants.py" "apps/SuffixTree.EmProto/$F" "$HERE/e5_${N}_mutants.json" "$HERE/results_e5_${TAG}${N}.tsv" "${SL:-}"
done
echo audit-done > "$HERE/e5_${TAG}audit.done"
