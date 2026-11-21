#!/bin/bash
# Inputs of work package 1b.9d whose LCP is wider than one byte or has a long overflow list (the growth of the container behind its reserved length is what the package reserves): the random-DNA files of the
# benchmark, with a part of the text repeated at its end. Every file has exactly as many symbols as its name says (one byte per symbol).
#   u16_<n>: the first 90 % of the text, then its first 10 %  -> LCP of two bytes per entry with an overflow list of about n / 10 pairs (n = 10^7: 9.9 M pairs)
#   u32_<n>: the first half of the text twice                  -> LCP of four bytes per entry (more than n / 8 values at or above 65,535)
# usage: make_1b9d_inputs.sh <dna text 10^7> <dna text 10^8> <output directory>
set -eu
D10=$1; D100=$2; OUT=$3
mkdir -p "$OUT"
for spec in "10m:$D10:10000000" "100m:$D100:100000000"; do
  IFS=: read -r name src n <<< "$spec"
  nine=$((n / 10 * 9)); tenth=$((n / 10)); half=$((n / 2))
  { head -c $nine "$src"; head -c $tenth "$src"; } > "$OUT/u16_$name.txt"
  { head -c $half "$src"; head -c $half "$src"; } > "$OUT/u32_$name.txt"
  for f in u16 u32; do [ "$(stat -c %s "$OUT/${f}_$name.txt")" = "$n" ] || { echo "wrong size of $f $name" >&2; exit 1; }; done
done
