#!/usr/bin/env bash
# ADR-005 section 5: the synthetic width measurement (chase / thr / scan / null) for w4, w5, w8, hyb, p34 at n = 10^7 .. 10^9, one and four threads, 5 repetitions.
# usage: run_width_synthetic.sh [sizes comma list]   Run alone on a quiet machine (the script prints the load); about 25-40 minutes for the full grid.
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
SIZES=${1:-10000000,100000000,500000000,1000000000}
R="$S2/results/width_synthetic.tsv"
"$S2/scripts/drift.sh" synthetic-start
quiet
$WB selftest | tail -1
$WB null --layouts w4,w5,w8,hyb,p34 --reps 7 --ops 20000000 --out "$R" --tag "null"
for n in ${SIZES//,/ }; do
  quiet
  $WB synthetic --layouts w4,w5,w8,hyb,p34 --n "$n" --v 6200000000 --threads 1,4 --reps 5 --ops 6000000 --out "$R" --tag "anon n=$n V=6.2e9"
done
"$S2/scripts/drift.sh" synthetic-end
