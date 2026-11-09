#!/usr/bin/env bash
# Work package 1b.4: the whole measurement campaign (nothing else may run on the machine). Raw rows in nav_raw.tsv (NAV_RAW overrides).
# MEM / MUM / anchors on dna_10m and rep_10m (5 fresh processes each, Ukkonen tree in the same process), the LCS family on two 1 Mbp texts, ns per link, the MEM time split, dna_50m MEM (information only).
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
D=${1:-/tmp/stbench}
RUNS=${2:-5}
export DOTNET_ROOT=${DOTNET_ROOT:-$HOME/.dotnet} PATH=$HOME/.dotnet:$PATH
dotnet build-server shutdown >/dev/null 2>&1
for X in dna_10m rep_10m; do
  for OP in mem mum anchors; do "$B/run_nav.sh" "$D" "$X" "$OP" "$RUNS"; done
done
"$B/run_nav.sh" "$D" dna_1m lcs "$RUNS" --other-kind mutated
RUN_BASE=100 "$B/run_nav.sh" "$D" dna_1m lcs "$RUNS" --other-kind random
for X in dna_10m rep_10m; do
  "$B/run_nav.sh" "$D" "$X" links "$RUNS"
  "$B/run_nav.sh" "$D" "$X" split 3
done
"$B/run_nav.sh" "$D" dna_50m mem 3 --ukkonen false
echo "# campaign finished $(date -u +%FT%TZ)" >> "${NAV_RAW:-$B/nav_raw.tsv}"
