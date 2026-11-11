#!/bin/bash
# Round 2 of the E3b audit: the non-equivalent survivors of round 1 and the two mutants that did not compile, against the tests added after round 1 (same runner, same watchdog).
#   usage: run_e3b_audit_round2.sh <worktree at the commit that contains the added tests>
set -u
W=${1:-/tmp/e3b-audit}
HERE=$(cd "$(dirname "$0")" && pwd)
export E3B_ROOT=$W PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1
for item in Dcs/SampleNamer.cs:samplenamer Dcs/ExternalDcs.cs:externaldcs Experiments/E3b.cs:e3b; do
  F=${item%%:*}; N=${item##*:}
  python3 -I "$HERE/run_e3b_mutants.py" "apps/SuffixTree.EmProto/$F" "$HERE/e3b_round2_${N}_mutants.json" "$HERE/results_e3b_round2_${N}.tsv"
done
echo audit-done > "$HERE/e3b_audit_round2.done"
