#!/usr/bin/env bash
# ADR-005 section 5: the dependent-read (chase) rows. The first version of the chain followed a random function and fell into its cycle of about sqrt(n) entries (cache-resident), so those
# rows were removed from results/width_synthetic.tsv; this version restarts the chain every 1,024 steps. usage: run_width_chase.sh [sizes]; run alone on a quiet machine.
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
SIZES=${1:-10000000,100000000,500000000,1000000000}
R="$S2/results/width_synthetic.tsv"
"$S2/scripts/drift.sh" chase-start
for n in ${SIZES//,/ }; do
  quiet
  $WB synthetic --layouts w4,w5,w8,hyb,p34 --n "$n" --v 6200000000 --threads 1,4 --kinds chase --reps 7 --ops 3000000 --out "$R" --tag "anon n=$n V=6.2e9 chain restarts every 1024"
done
"$S2/scripts/drift.sh" chase-end
