#!/bin/bash
# E7 setup (outside the shipped code): for every input name <family>_<size>m (dna_50m, dna_100m, dna_200m, dna_500m, genome_50m, ...) the text file (the input of the pipeline), the references of gate S3 made by the in-RAM
# oracles of the earlier stages (as run_e6_setup.sh), and the SHA-256 of every reference file:
#   <name>.txt            ASCII text, one byte per symbol (e6-mktext)
#   <name>.ref.pk         packed text of the in-RAM `pack` command           <name>.refsa           suffix array of SaisKernel in RAM (oracle --out)
#   <name>.oracle.lcp/.ovf compact LCP of LcpKernel on that SA (lcp-oracle)   <name>.ref.cldP/.stats  child table + statistics of the in-RAM ChildTable.Build / EsaStatistics (e5-oracle)
#   <name>.refs.sha       "<kind> <sha256>" per reference file (pk sa lcp ovf cld) and the four statistics the pipeline compares ("stat key=value"); the timed runs compare their outputs with these hashes
# 100 Mbp inputs: the files of E6 (E7_LINK_FROM, default /tmp/stbench/e6) are linked, not regenerated.
# The disk is tight at 500 Mbp: E7_DROP_SA=1 deletes the suffix array reference as soon as its hash and the LCP files are made (e5-oracle does not read it), E7_DROP_REFS=1 deletes all the big reference files at the end
# (the text and the hashes stay). Every step is marked by <name>.step.<kind> when it is complete, so an interrupted step is redone from scratch.
# usage: run_e7_setup.sh <work dir> <rows file> <names...>
set -u
W=$1; ROWS=$2; shift 2
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e7_run/SuffixTree.EmProto.dll}
DATA=${E7_DATA:-/tmp/stbench/data}
LINK=${E7_LINK_FROM:-/tmp/stbench/e6}
cd "$W" || exit 1
printf 'R\tE7\tsetup\thead\t%s\t\tgit rev-parse HEAD of the checkout when the inputs were made\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" >> "$ROWS"

hash1() { sha256sum "$1" | cut -d' ' -f1; }

for X in "$@"; do
  FAM=${X%%_*}; SZ=${X#*_}; N=$(( ${SZ%m} * 1000000 )); F="$FAM:$N"
  unset DOTNET_GCHeapHardLimit
  P=$X.refs.part
  if [ ! -f "$X.refs.sha" ]; then
    if [ -f "$LINK/$X.ref.cldP" ] && [ -f "$LINK/$X.txt" ]; then
      for S in txt ref.pk ref.pk.meta refsa oracle.lcp oracle.ovf ref.cldP ref.stats; do [ -e "$X.$S" ] || ln -s "$LINK/$X.$S" "$X.$S"; done
      touch "$X.step.txt" "$X.step.pk" "$X.step.sa" "$X.step.lcp" "$X.step.cld"
    fi
    [ -f "$X.step.txt" ] || { rm -f "$X.txt"; dotnet "$DLL" e6-mktext --family "$F" --data "$DATA" --out "$X.txt" --case "$X" --rows "$ROWS" > /dev/null && touch "$X.step.txt"; }
    [ -f "$X.step.pk" ] || { rm -f "$X.ref.pk"; dotnet "$DLL" pack --family "$F" --data "$DATA" --out "$X.ref.pk" --rows "$ROWS" > /dev/null && touch "$X.step.pk"; }
    if [ ! -f "$X.step.sa" ]; then
      rm -f "$X.refsa"; dotnet "$DLL" oracle --family "$F" --data "$DATA" --out "$X.refsa" --rows "$ROWS" > /dev/null && touch "$X.step.sa"
    fi
    if [ ! -f "$X.step.lcp" ]; then
      [ -f "$X.refsa" ] || { echo "$X: the suffix array reference is gone before its LCP exists; remove $X.step.sa and rerun" >&2; exit 3; }
      rm -f "$X.oracle.lcp" "$X.oracle.ovf"; dotnet "$DLL" lcp-oracle --text "$X.ref.pk" --sa "$X.refsa" --out "$X.oracle" --case "$X" --rows "$ROWS" > /dev/null && touch "$X.step.lcp"
    fi
    grep -q '^sa ' "$P" 2>/dev/null || { [ -f "$X.refsa" ] && echo "sa $(hash1 "$X.refsa")" >> "$P"; }
    if [ "${E7_DROP_SA:-0}" = 1 ]; then rm -f "$X.refsa"; fi
    if [ ! -f "$X.step.cld" ]; then
      rm -f "$X.ref.cldP" "$X.ref.cldT" "$X.ref.psv" "$X.ref.nsv" "$X.ref.nsym" "$X.ref.stats"
      dotnet "$DLL" e5-oracle --text "$X.ref.pk" --lcp "$X.oracle" --out "$X.ref" --case "$X" --rows "$ROWS" > /dev/null && touch "$X.step.cld"
      rm -f "$X.ref.cldT" "$X.ref.psv" "$X.ref.nsv" "$X.ref.nsym"
    fi
    for K in "pk $X.ref.pk" "lcp $X.oracle.lcp" "ovf $X.oracle.ovf" "cld $X.ref.cldP"; do set -- $K; grep -q "^$1 " "$P" 2>/dev/null || echo "$1 $(hash1 "$2")" >> "$P"; done
    grep -q '^sa ' "$P" || { echo "$X: no hash of the suffix array reference" >&2; exit 3; }
    {
      for K in pk sa lcp ovf cld; do grep "^$K " "$P"; done
      grep -E '^(nodecount|internal|lrsdepth|lrsrank)=' "$X.ref.stats" | sort | sed 's/^/stat /'
    } > "$X.refs.sha.tmp" && mv "$X.refs.sha.tmp" "$X.refs.sha"
    printf 'R\tE7\t%s\tsetup\treference_sizes\t%s\tbytes\tpk lcp ovf cld of the references (the suffix array has 4 (n + 1) bytes)\n' "$X" "$(stat -L -c %s "$X.ref.pk" "$X.oracle.lcp" "$X.oracle.ovf" "$X.ref.cldP" 2>/dev/null | tr '\n' ' ')" >> "$ROWS"
  fi
  if [ "${E7_DROP_REFS:-0}" = 1 ]; then rm -f "$X.ref.pk" "$X.refsa" "$X.oracle.lcp" "$X.oracle.ovf" "$X.ref.cldP"; fi
  echo "$X ready"
done
echo setup-done
