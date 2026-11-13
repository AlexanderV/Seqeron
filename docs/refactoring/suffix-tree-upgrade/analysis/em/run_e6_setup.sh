#!/bin/bash
# E6 setup (outside the shipped code): for every input name the text file (the input of the pipeline), the references of gate G1 made by the earlier stages' in-RAM oracles:
#   <name>.ref.pk         packed text of the in-RAM `pack` command (the reference of PH0)
#   <name>.refsa          suffix array of SaisKernel in RAM (oracle --out): n + 1 little-endian int32 (the reference of PH1-PH3)
#   <name>.oracle.lcp/.ovf  compact LCP of LcpKernel on that suffix array (lcp-oracle; the reference of PH4)
#   <name>.ref.cldP, .ref.stats  child table of the packed layout of the product's ChildTable.Build and the statistics (e5-oracle; the reference of PH5)
# usage: run_e6_setup.sh <work dir> <rows file> <names...>   names: see family_of below
set -u
W=$1; ROWS=$2; shift 2
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e6_run/SuffixTree.EmProto.dll}
DATA=${E6_DATA:-/tmp/stbench/data}
cd "$W" || exit 1
printf 'R\tE6\tsetup\thead\t%s\t\tgit rev-parse HEAD of the checkout when the inputs were made\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" >> "$ROWS"

family_of() { # name -> family spec (a file input of /tmp/stbench/data has its own name)
  case $1 in
    dna_100m) echo dna:100000000 ;;
    genome_100m) echo genome:100000000 ;;
    dna_10m|rep_10m|prot_1m|dna_1m) echo "$1" ;;
    a50) echo a:50000000 ;;
    nrun30) echo nrun:30000000:20000000 ;;
    fib10) echo fib:10000000 ;;
    fam_a_1e5) echo a:100000 ;;
    fam_a_1e6) echo a:1000000 ;;
    fam_fib_1e5) echo fib:100000 ;;
    fam_fib_1e6) echo fib:1000000 ;;
    fam_tm_1e5) echo tm:100000 ;;
    fam_tm_1e6) echo tm:1000000 ;;
    fam_tandem_1e5) echo tandem:171x600 ;;
    fam_tandem_1e6) echo tandem:171x6000 ;;
    fam_nrun_1e5) echo nrun:100000:10000 ;;
    fam_nrun_1e6) echo nrun:1000000:100000 ;;
    fam_genome_1e6) echo genome:1000000 ;;
    *) echo "unknown $1" >&2; return 1 ;;
  esac
}

for X in "$@"; do
  F=$(family_of "$X") || exit 2
  unset DOTNET_GCHeapHardLimit
  case $F in
    dna_10m|rep_10m|prot_1m|dna_1m) [ -e "$X.txt" ] || ln -s "$DATA/$F.txt" "$X.txt" ;;
    *) [ -f "$X.txt" ] || dotnet "$DLL" e6-mktext --family "$F" --data "$DATA" --out "$X.txt" --case "$X" --rows "$ROWS" > /dev/null ;;
  esac
  [ -f "$X.ref.pk" ] || dotnet "$DLL" pack --family "$F" --data "$DATA" --out "$X.ref.pk" --rows "$ROWS" > /dev/null
  [ -f "$X.refsa" ] || dotnet "$DLL" oracle --family "$F" --data "$DATA" --out "$X.refsa" --rows "$ROWS" > /dev/null
  [ -f "$X.oracle.lcp" ] || dotnet "$DLL" lcp-oracle --text "$X.ref.pk" --sa "$X.refsa" --out "$X.oracle" --case "$X" --rows "$ROWS" > /dev/null
  if [ ! -f "$X.ref.cldP" ]; then
    dotnet "$DLL" e5-oracle --text "$X.ref.pk" --lcp "$X.oracle" --out "$X.ref" --case "$X" --rows "$ROWS" > /dev/null
    rm -f "$X.ref.cldT" "$X.ref.psv" "$X.ref.nsv" "$X.ref.nsym"
  fi
  echo "$X ready"
done
echo setup-done
