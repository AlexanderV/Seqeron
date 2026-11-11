#!/bin/bash
# Round 3 of the E4 audit: the three mutants that survived round 2 (B06, B06b, E08b) against the two tests that were strengthened after it (same runner, same watchdog).
#   usage: run_e4_audit_round3.sh <worktree at the commit that contains the strengthened tests>
set -u
W=${1:-/tmp/e4-audit}
HERE=$(cd "$(dirname "$0")" && pwd)
export E4_ROOT=$W PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
for item in Ph4/LcpBuild.cs:build Experiments/E4.cs:e4; do
  F=${item%%:*}; N=${item##*:}
  python3 -I "$HERE/run_e4_mutants.py" "apps/SuffixTree.EmProto/$F" "$HERE/e4_round3_${N}_mutants.json" "$HERE/results_e4_round3_${N}.tsv"
done
echo audit-done > "$HERE/e4_audit_round3.done"
