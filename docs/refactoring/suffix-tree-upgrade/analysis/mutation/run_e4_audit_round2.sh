#!/bin/bash
# Round 2 of the E4 audit: the non-equivalent survivors of round 1, the replacements of the mutants that did not compile and the mutants of the exit-code helper, against the tests added after round 1
# (same runner, same watchdog).   usage: run_e4_audit_round2.sh <worktree at the commit that contains the added tests>
set -u
W=${1:-/tmp/e4-audit}
HERE=$(cd "$(dirname "$0")" && pwd)
export E4_ROOT=$W PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
for item in Ph4/LcpBuild.cs:build Ph4/PlcpScan.cs:plcpscan Experiments/E4.cs:e4; do
  F=${item%%:*}; N=${item##*:}
  python3 -I "$HERE/run_e4_mutants.py" "apps/SuffixTree.EmProto/$F" "$HERE/e4_round2_${N}_mutants.json" "$HERE/results_e4_round2_${N}.tsv"
done
echo audit-done > "$HERE/e4_audit_round2.done"
