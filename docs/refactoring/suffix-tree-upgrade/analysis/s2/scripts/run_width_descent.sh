#!/usr/bin/env bash
# ADR-005 section 5: the count / locate descent over a real container with the stored width swapped (widthbench `descent`), the variants interleaved.
# usage: run_width_descent.sh <container> <label> <variants> [reps] [queries]
#   e.g.  run_width_descent.sh $S2_OUT/work_dna100/out.sesa dna_100m v0,v1,v1b,v2,v3,v4,v5,v6,v7 11 200000
# The container is built by `SuffixTree.Esa.Bench sesa build` (README.md of analysis/s2). Run alone on a quiet machine; the first variant is checked against the product reader (--check).
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
C=${1:?container}; L=${2:?label}; V=${3:?variants}; REPS=${4:-11}; Q=${5:-200000}
R="$S2/results/width_descent.tsv"
quiet
"$S2/scripts/drift.sh" "descent-$L-start"
$WB descent --container "$C" --variants "$V" --threads 1,4 --reps "$REPS" --queries "$Q" --check true --out "$R" --tag "$L"
"$S2/scripts/drift.sh" "descent-$L-end"
