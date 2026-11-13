#!/bin/bash
# E5 setup (outside the shipped code): the packed texts, the in-RAM rank files, the suffix array files of the E3 builder and the LcpKernel files (the compact LCP in the two-file form of PH4: <name>.oracle.lcp / .oracle.ovf,
# which is the input form of PH5), for the inputs of the E5 protocol. Same method as run_e4_setup.sh.
# usage: run_e5_setup.sh <work dir> <rows file> <names...>   names: see the case list below; E3 / E4 files of the earlier prototypes (dna_100m, genome_100m) are taken from /tmp/stbench/e3 and /tmp/stbench/orch4 (read only)
# Every SA build runs in a fresh process after sync + drop_caches with DOTNET_GCHeapHardLimit = 96 MiB, exactly as the E3 protocol.
set -u
W=$1; ROWS=$2; shift 2
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e5_run/SuffixTree.EmProto.dll}
cd "$W" || exit 1
mkdir -p scratch
printf 'R\tE5\tsetup\thead\t%s\t\tgit rev-parse HEAD of the checkout when the inputs were made\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" >> "$ROWS"

family_of() { # name -> family spec for the pack command, or "-" when the packed file exists
  case $1 in
    a10) echo a:10000000 ;;
    a25) echo a:25000000 ;;
    a50) echo a:50000000 ;;
    a200) echo a:200000000 ;;
    fib10) echo fib:10000000 ;;
    fib30) echo fib:30000000 ;;
    nrun30) echo nrun:30000000:20000000 ;;
    dna10) echo dna_10m ;;
    rep10) echo rep_10m ;;
    dna10r) echo dna:10000000 ;;
    dna50r) echo dna:50000000 ;;
    genome10) echo genome:10000000 ;;
    tandem) echo tandem:171x10000 ;;
    nrun10) echo nrun:10000000:6000000 ;;
    dna_100m|genome_100m) echo - ;;
    *) echo "unknown $1" >&2; return 1 ;;
  esac
}

for X in "$@"; do
  F=$(family_of "$X") || exit 2
  if [ "$F" = "-" ]; then
    for ext in pk pk.meta ranks ranks.meta; do [ -e "$X.$ext" ] || ln -s "/tmp/stbench/e3/$X.$ext" "$X.$ext"; done
    if [ "$X" = dna_100m ]; then
      for ext in sa oracle.lcp oracle.ovf; do [ -e "$X.$ext" ] || cp "/tmp/stbench/orch4/$X.$ext" "$X.$ext"; done  # copies, not links: the programs take the length of the file from FileInfo, which is the length of a link
    fi
  else
    ( unset DOTNET_GCHeapHardLimit
      [ -f "$X.pk" ] || dotnet "$DLL" pack --family "$F" --data /tmp/stbench/data --out "$X.pk" > /dev/null
      [ -f "$X.ranks" ] || dotnet "$DLL" mkranks --text "$X.pk" --out "$X.ranks" > /dev/null )
  fi
  export DOTNET_GCHeapHardLimit=0x6000000
  if [ ! -e "$X.sa" ]; then
    sync; echo 3 > /proc/sys/vm/drop_caches
    dotnet "$DLL" e3 --text "$X.pk" --ranks "$X.ranks" --sa "$X.sa" --scratch scratch --threads 4 --case "$X.sa_phase" --rows "$ROWS" > /dev/null
    echo "$X SA rc=$?"
  fi
  unset DOTNET_GCHeapHardLimit
  [ -f "$X.oracle.lcp" ] || dotnet "$DLL" lcp-oracle --text "$X.pk" --sa "$X.sa" --out "$X.oracle" --case "$X" --rows "$ROWS" > /dev/null
  echo "$X ready"
done
echo setup-done
