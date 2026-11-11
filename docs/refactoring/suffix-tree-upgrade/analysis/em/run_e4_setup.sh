#!/bin/bash
# E4 setup (outside the shipped code): the packed texts, the in-RAM rank files, the suffix array files of the E3 builder (the SA phase is timed here, cold cache, 4 threads, the reference of the abandon rule
# "LCP phase > 3x the SA phase") and the LcpKernel oracle files, for the inputs of the E4 protocol.
# usage: run_e4_setup.sh <work dir> <rows file> <names...>   names: see the case list below; E3 files of the earlier prototypes (dna_100m, genome_100m) are taken from /tmp/stbench/e3
# Every SA build runs in a fresh process after sync + drop_caches with DOTNET_GCHeapHardLimit = 96 MiB, exactly as the E3 protocol.
set -u
W=$1; ROWS=$2; shift 2
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e4_run/SuffixTree.EmProto.dll}
SA_RUNS=${SA_RUNS:-1}
cd "$W" || exit 1
mkdir -p scratch
printf 'R\tE4\tsetup\thead\t%s\t\tgit rev-parse HEAD of the binary under test\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" >> "$ROWS"

family_of() { # name -> family spec for the pack command, or "-" when the packed file exists
  case $1 in
    dna10) echo dna_10m ;;
    rep10) echo rep_10m ;;
    a50) echo a:50000000 ;;
    a10) echo a:10000000 ;;
    fib10) echo fib:10000000 ;;
    tandem) echo tandem:171x10000 ;;
    nrun10) echo nrun:10000000:1000000 ;;
    genome10) echo genome:10000000 ;;
    genome50) echo genome:50000000 ;;
    dna10r) echo dna:10000000 ;;
    dna50r) echo dna:50000000 ;;
    dna200r) echo dna:200000000 ;;
    prot1m) echo prot_1m ;;
    dna1m) echo dna_1m ;;
    fam_a_1e5) echo a:100000 ;;
    fam_a_1e6) echo a:1000000 ;;
    fam_fib_1e5) echo fib:100000 ;;
    fam_fib_1e6) echo fib:1000000 ;;
    fam_tm_1e5) echo tm:100000 ;;
    fam_tm_1e6) echo tm:1000000 ;;
    fam_tandem_1e5) echo tandem:171x600 ;;
    fam_nrun_1e6) echo nrun:1000000:100000 ;;
    dna_100m|genome_100m) echo - ;;
    *) echo "unknown $1" >&2; return 1 ;;
  esac
}

for X in "$@"; do
  F=$(family_of "$X") || exit 2
  if [ "$F" = "-" ]; then
    for ext in pk pk.meta ranks ranks.meta; do [ -e "$X.$ext" ] || ln -s "/tmp/stbench/e3/$X.$ext" "$X.$ext"; done
  else
    ( unset DOTNET_GCHeapHardLimit
      [ -f "$X.pk" ] || dotnet "$DLL" pack --family "$F" --data /tmp/stbench/data --out "$X.pk" > /dev/null
      [ -f "$X.ranks" ] || dotnet "$DLL" mkranks --text "$X.pk" --out "$X.ranks" > /dev/null )
  fi
  export DOTNET_GCHeapHardLimit=0x6000000
  for i in $(seq 1 "$SA_RUNS"); do
    sync; echo 3 > /proc/sys/vm/drop_caches
    dotnet "$DLL" e3 --text "$X.pk" --ranks "$X.ranks" --sa "$X.sa" --scratch scratch --threads 4 --case "$X.sa_phase" --run "$i" --rows "$ROWS" > /dev/null
    echo "$X SA run $i rc=$?"
  done
  unset DOTNET_GCHeapHardLimit
  [ -f "$X.oracle.lcp" ] || dotnet "$DLL" lcp-oracle --text "$X.pk" --sa "$X.sa" --out "$X.oracle" --case "$X" --rows "$ROWS" > /dev/null
  echo "$X ready"
done
echo setup-done
