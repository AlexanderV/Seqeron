#!/bin/bash
# E5 protocol (outside the shipped code): PH5 (the lcp-interval scan with the external stack, layouts A and F) on the inputs prepared by run_e5_setup.sh.
# usage: run_e5_protocol.sh <work dir with <name>.pk, <name>.sa, <name>.oracle.lcp/.ovf> <rows file> <group> [names...]
# groups: adv   (G1/G3 on a^n, N run, Fibonacci, the 100 Mbp inputs: layouts a, f, both; 3 runs each, the first run compared byte for byte with the in-RAM reference, the others not)
#         shm   (the same on a^n / N run / Fibonacci with the spill file on tmpfs: the isolation run of the write amplification)
#         slope (RssAnon over n: layout both, 3 runs)
#         fused (PH4 + PH5 in one process through the B3 hook, scale-down emulation of PH4 as in E4; then the NSV pass)
#         ph4ref (the E4 base mode on the same input, the reference of "PH5 time must not exceed the PH4 time")
#         nosp  (spill file not flushed: the dirty pages of the spill, information)
#         par   (layout f and both with the NSV pass on its own thread while the forward scan runs: --parallel-nsv)
#         a200  (a^n 200 M, layout a only, no in-RAM reference: the output of the last run is compared with the closed form by check_child_allsame.py; inputs from run_e5_setup.sh a200)
#         stats (the scan alone: no output file, the spill on tmpfs: the cost of the stack and the scan without any write behind the cursor)
#         slots (layout a with 1 / 2 / 4 windows instead of 8: how close to thrashing the writes behind the cursor are)
#         nextsym (the optional NEXTSYM pass after the statistics-only scan, report only; the reference .nsym from the in-RAM oracle with --sa)
# Every run: dropped page cache (cold), fresh process, DOTNET_GCHeapHardLimit = 96 MiB, real block sizes (256 K entries, 1 MiB windows), the independent sampler sample_proc.py, the block device's sectors written before the
# run and after the run + sync (definition 2 of the write amplification).
set -u
W=$1; ROWS=$2; GROUP=$3; shift 3
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e5_run/SuffixTree.EmProto.dll}
HERE=$(cd "$(dirname "$0")" && pwd)
DEV=${E5_DEV:-vda}
cd "$W" || exit 1
mkdir -p scratch /dev/shm/e5spill
printf 'R\tE5\tprotocol\thead\t%s\t\tgit rev-parse HEAD of the binary under test; group %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$GROUP" >> "$ROWS"

sectors() { awk '{print $7}' /sys/block/$DEV/stat; }

run() { # input label run mode(file|fused|ph4) extra-args...
  T=$1; L=$2; R=$3; MODE=$4; shift 4
  rm -f $T.e5.* /dev/shm/e5spill/* scratch/* 2>/dev/null
  rm -f $T.pidfile
  sync; echo 3 > /proc/sys/vm/drop_caches
  python3 -I "$HERE/sample_proc.py" $T.pidfile "$ROWS" "$L" "$R" &
  SP=$!
  sleep 0.2
  S0=$(sectors)
  case $MODE in
    file)  ( echo $BASHPID > $T.pidfile; exec dotnet "$DLL" e5 --text $T.pk --lcp $T.oracle --out $T.e5 --scratch scratch --case "$L" --run "$R" --rows "$ROWS" "$@" > /dev/null ) ;;
    fused) ( echo $BASHPID > $T.pidfile; exec dotnet "$DLL" e5 --fused --text $T.pk --sa $T.sa --out $T.e5 --scratch scratch --scale-n 3000000000 --case "$L" --run "$R" --rows "$ROWS" "$@" > /dev/null ) ;;
    ph4)   ( echo $BASHPID > $T.pidfile; exec dotnet "$DLL" e4 --text $T.pk --sa $T.sa --out $T.p4 --scratch scratch --scale-n 3000000000 --case "$L" --run "$R" --rows "$ROWS" "$@" > /dev/null ) ;;
  esac
  RC=$?
  wait $SP
  sync
  S1=$(sectors)
  printf 'R\tE5\t%s\t%s\tdevice_write_bytes\t%s\tbytes\tsectors written by %s during the run and the sync that follows it (512 B each; includes the spill file, journal, everything)\n' "$L" "$R" "$(( (S1 - S0) * 512 ))" "$DEV" >> "$ROWS"
  echo "$L run $R rc=$RC"
  [ -n "${E5_KEEP:-}" ] || rm -f $T.e5.* $T.p4.* /dev/shm/e5spill/* 2>/dev/null
}

oracle() { # input: the in-RAM reference files for --verify
  T=$1
  ( unset DOTNET_GCHeapHardLimit; dotnet "$DLL" e5-oracle --text $T.pk --lcp $T.oracle --out $T.ref --case "$T.oracle" --rows "$ROWS" > /dev/null )
}

oracle_sa() { # input: the in-RAM reference files including NEXTSYM
  T=$1
  ( unset DOTNET_GCHeapHardLimit; dotnet "$DLL" e5-oracle --text $T.pk --lcp $T.oracle --sa $T.sa --out $T.ref --case "$T.oracle" --rows "$ROWS" > /dev/null )
}

set_names() { NAMES=("$@"); }
case $GROUP in
adv)
  [ $# -gt 0 ] && set_names "$@" || set_names a10 a50 nrun30 fib10 fib30 rep10 dna10 dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    oracle $T
    for LAY in a f both; do
      run $T "$T.$LAY" 1 file --layout $LAY --verify $T.ref
      for i in 2 3; do run $T "$T.$LAY" $i file --layout $LAY; done
    done
    run $T "$T.apacked" 1 file --layout a --codec packed --verify $T.ref
    rm -f $T.ref.*
  done ;;
shm)
  [ $# -gt 0 ] && set_names "$@" || set_names a50 nrun30 fib30 dna_100m
  for T in "${NAMES[@]}"; do
    for LAY in a f; do
      for i in 1 2 3; do run $T "$T.$LAY.shm" $i file --layout $LAY --spill-dir /dev/shm/e5spill; done
    done
  done ;;
nosp)
  [ $# -gt 0 ] && set_names "$@" || set_names a50 nrun30
  for T in "${NAMES[@]}"; do
    run $T "$T.a.nospillflush" 1 file --layout a --spill-flush-mb 0
    run $T "$T.a.noflush" 1 file --layout a --spill-flush-mb 0 --no-flush
  done ;;
slope)
  [ $# -gt 0 ] && set_names "$@" || set_names a10 a25 a50 dna10r dna50r dna_100m
  for T in "${NAMES[@]}"; do
    for LAY in both; do
      for i in 1 2 3; do run $T "$T.$LAY.slope" $i file --layout $LAY; done
      for i in 1 2 3; do run $T "$T.$LAY.slopetouch" $i file --layout $LAY --touch-buffers; done
    done
  done ;;
fused)
  [ $# -gt 0 ] && set_names "$@" || set_names dna_100m genome_100m a50
  for T in "${NAMES[@]}"; do
    oracle $T
    run $T "$T.fusedA" 1 fused --layout a --verify $T.ref
    run $T "$T.fusedF" 1 fused --layout f --verify $T.ref
    for i in 2 3; do run $T "$T.fusedA" $i fused --layout a; run $T "$T.fusedF" $i fused --layout f; done
    rm -f $T.ref.*
  done ;;
stats)
  [ $# -gt 0 ] && set_names "$@" || set_names a50 nrun30 fib30 dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    for i in 1 2 3; do run $T "$T.stats.shm" $i file --layout stats --spill-dir /dev/shm/e5spill; done
  done ;;
slots)
  [ $# -gt 0 ] && set_names "$@" || set_names a50 nrun30 fib30 genome_100m dna_100m
  for T in "${NAMES[@]}"; do
    for SL in 1 2 4; do run $T "$T.a.slots$SL" 1 file --layout a --slots $SL; done
  done ;;
par)
  [ $# -gt 0 ] && set_names "$@" || set_names a50 nrun30 fib30 dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    for LAY in f both; do
      for i in 1 2 3; do run $T "$T.$LAY.par" $i file --layout $LAY --parallel-nsv; done
    done
  done ;;
a200)
  T=a200
  N=$(grep '^n=' $T.pk.meta | cut -d= -f2)
  for i in 1 2; do run $T "$T.a" $i file --layout a; done
  E5_KEEP=1 run $T "$T.a" 3 file --layout a
  python3 -I "$HERE/check_child_allsame.py" $T.e5.cld $N | sed "s/^/a200 closed form: /"
  rm -f $T.e5.* ;;
nextsym)
  [ $# -gt 0 ] && set_names "$@" || set_names dna10 a50 nrun30 dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    oracle_sa $T
    run $T "$T.nextsym" 1 file --layout stats --nextsym --sa $T.sa --verify $T.ref
    for i in 2 3; do run $T "$T.nextsym" $i file --layout stats --nextsym --sa $T.sa; done
    rm -f $T.ref.*
  done ;;
ph4ref)
  [ $# -gt 0 ] && set_names "$@" || set_names a10 a50 nrun30 fib10 fib30 rep10 dna10 dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    for i in 1 2 3; do run $T "$T.ph4" $i ph4 --threads 4; done
  done ;;
*) echo "unknown group $GROUP"; exit 2 ;;
esac
echo protocol-done
