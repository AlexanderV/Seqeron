#!/usr/bin/env bash
# ADR-005 section 5: the descent at n = 5x10^8 (the page-cache and DRAM regime of this machine): the variants in groups that fit in 15 GB of RAM, each group interleaved; v1 is the anchor of every group
# that has room for it (v5, the 8-byte layout, needs 8 GB alone). usage: run_width_descent_500m.sh <container> [reps] [queries]   (the container of a 500 Mbp text, Core profile is enough)
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
C=${1:?container}; REPS=${2:-9}; Q=${3:-200000}
R="$S2/results/width_descent.tsv"
"$S2/scripts/drift.sh" descent-500m-start
for g in "v0,v1,v1b,v2,v8h,v9h" "v1,v3" "v1,v6" "v1,v7" "v5"; do
  quiet
  echo "== group $g"
  $WB descent --container "$C" --variants "$g" --threads 1,4 --reps "$REPS" --queries "$Q" --check false --out "$R" --tag "dna_500m group $g"
done
"$S2/scripts/drift.sh" descent-500m-end
