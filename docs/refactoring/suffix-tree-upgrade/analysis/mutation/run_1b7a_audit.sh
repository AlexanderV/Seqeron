#!/bin/bash
# Runs the hand-made mutants of the production code of 1b.7a, one source file after the other, in a separate git worktree (AUDIT_ROOT) so that the working checkout is never mutated. Every build and test run has a
# watchdog (run_1b7a_mutants.py kills the whole process group; a hang counts as killed). Results: results_1b7a_<round><file>.tsv next to this script.
#   usage: run_1b7a_audit.sh <worktree> [round prefix: "", "r2_", "r3_"]
set -u
W=${1:-/tmp/1b7a-audit}
R=${2:-}
HERE=$(cd "$(dirname "$0")" && pwd)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
mkdir -p "$W/tests/SuffixTree/SuffixTree.Esa.Slim"
cp "$HERE/SuffixTree.Esa.Slim.csproj.template" "$W/tests/SuffixTree/SuffixTree.Esa.Slim/SuffixTree.Esa.Slim.csproj"
for L in "$HERE"/1b7a_${R}*_listed.json; do
  [ -f "$L" ] || continue
  N=$(basename "$L" _listed.json)
  [ -n "$R" ] || case "$N" in 1b7a_r2_*|1b7a_r3_*) continue ;; esac
  python3 -I "$HERE/run_1b7a_mutants.py" "$W" "$L" "$HERE/results_$N.tsv"
done
echo "audit-done $R" > "$HERE/1b7a_${R}audit.done"
