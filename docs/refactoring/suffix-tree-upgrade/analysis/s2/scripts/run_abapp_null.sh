#!/usr/bin/env bash
# The A/B procedure that the gates of S2.3-S2.5 use (ADR-005 section 3.4), run here as its null experiment: the same build against itself on a 10^8 container, so that the noise floor of the instrument is on record before any change.
# For a real A/B, specA and specB name two build directories (the S1 head built to its own directory, the S2 head to another): dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o <dir>.
# usage: run_abapp_null.sh <esa build dir> <container> <pattern-source dir with <name>.present20.txt, .absent20.txt, .locate8.txt> <name> [threads] [processes]
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
ESA=${1:?esa dir}; C=${2:?container}; Q=${3:?queries dir}; N=${4:?name}; TH=${5:-4}; PROCS=${6:-3}
ABAPP_BIN=${ABAPP_BIN:-$D/b9/bin/abapp}                         # built from analysis/bench/abapp (S1: run_1b9_gate.sh header)
P=$OUT/ab_patterns_$N.txt
{ sed 's/^/P /' "$Q/$N.present20.txt"; sed 's/^/L /' "$Q/$N.locate8.txt"; sed 's/^/A /' "$Q/$N.absent20.txt"; } > "$P"
R=$S2/results/abapp_null.tsv
quiet
"$S2/scripts/drift.sh" "abapp-null-$N-start"
ABAPP="$ABAPP_BIN" "$ROOT/docs/refactoring/suffix-tree-upgrade/analysis/bench/abapp/run_ab.sh" "$R" "$PROCS" "$P" 7 "$TH" 1.5 "$ESA|file:$C" "$ESA|file:$C" "null-$N-t$TH"
"$S2/scripts/drift.sh" "abapp-null-$N-end"
