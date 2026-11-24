#!/usr/bin/env bash
# ADR-005 section 6: the read cost of the packed DNA text (2 and 3 bits per symbol) against the byte text on the 10 Mbp inputs of the S1 gate (dna_10m: random; rep_10m: a repeated block, long common extensions).
# usage: run_text_bench.sh [reps]   Run alone on a quiet machine.
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
REPS=${1:-9}
quiet
"$S2/scripts/drift.sh" text-start
$WB text --text "$D/data/dna_10m.txt" --rep "$D/data/rep_10m.txt" --reps "$REPS" --ops 20000000 --out "$S2/results/text_bench.tsv" --tag "10 Mbp"
"$S2/scripts/drift.sh" text-end
