#!/bin/bash
# E9 protocol (outside the shipped code): the timed runs of the Psi pass and of the routes it is compared with, on the inputs prepared by run_e9_setup.sh.
# usage: run_e9_protocol.sh <work dir> <rows file> <group> <names...>
# groups: psi      PH6, Psi by sigma-way distribution of the BWT (e9-psi)                                  label <name>.e9.psi
#         psi256   the same with all of a byte alphabet in one pass (--symbols-per-pass 256: 256 windows of 64 KB)   label <name>.e9.psi256
#         bwt      the BWT pass over the SA file and the packed text (e9-bwt)                               label <name>.e9.bwt
#         isa      ISA by DIST (e9-isa, --scale-n 3000000000: R = 8,000,000 n / 3x10^9, P = 375; the chunk write buffer 64 KB x n / 10^8, at least 1 KB)   label <name>.e9.isa
#         psiisa   Psi = ISA[SA[i] + 1] from the ISA file (e9-psi-isa; needs <name>.ref.isa)              label <name>.e9.psiisa
#         base     the cost of starting the runtime: e9-psi on a BWT of 1,001 rows (the names are ignored)  label e9.base
# Every run: sync and dropped page cache (cold), fresh process, DOTNET_GCHeapHardLimit = 96 MiB (0x6000000), the independent sampler sample_e9.py (RssAnon at 20 ms, Dirty, /proc/<pid>/io, /sys/block/vda/stat before and
# after the exit and after sync), the program's own rows, then (outside the timed region) cmp of the output with the reference, the output is deleted. A run whose rows already end in run_done is skipped (restart safe);
# an interrupted run's rows are purged and the run is repeated cold. E9_RUNS runs per cell (default 3).
set -u
W=$1; ROWS=$2; GROUP=$3; shift 3
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e9_run/SuffixTree.EmProto.dll}
HERE=$(cd "$(dirname "$0")" && pwd)
TIMEOUT=${E9_TIMEOUT:-1200}
cd "$W" || exit 1
printf 'R\tE9\tprotocol\thead\t%s\t\tgit rev-parse HEAD of the checkout when the group started; group %s; binary %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$GROUP" "$DLL" >> "$ROWS"

purge() { awk -F'\t' -v l="$1" -v r="$2" '!($3==l && $4==r)' "$ROWS" > "$ROWS.tmp" && mv "$ROWS.tmp" "$ROWS"; }

run() { # input label run scratch-dir reference-file output-file command-args...
  T=$1; L=$2; R=$3; SCR=$4; REF=$5; OUT=$6; shift 6
  if grep -q -P "^R\tE9\t$L\t$R\trun_done\t" "$ROWS"; then echo "$L run $R already done"; return 0; fi
  purge "$L" "$R"
  PIDF=pid.$L.$R
  rm -rf "$SCR" "$OUT" $PIDF 2>/dev/null
  FREE=$(df -B1G --output=avail /tmp | tail -1 | tr -d ' ')
  printf 'R\tE9\t%s\t%s\tdisk_free_at_start\t%s\tGB\tdf of /tmp before the run\n' "$L" "$R" "$FREE" >> "$ROWS"
  sync; echo 3 > /proc/sys/vm/drop_caches
  python3 -I "$HERE/sample_e9.py" $PIDF "$ROWS" "$L" "$R" "$SCR" &
  SP=$!
  sleep 0.2
  S0=$(date +%s.%N)
  dotnet "$DLL" "$@" --case "$L" --run "$R" --rows "$ROWS.prog" > /dev/null &
  PID=$!
  echo $PID > $PIDF
  ( exec > /dev/null 2>&1; for ((w = 0; w < TIMEOUT; w++)); do kill -0 $PID 2>/dev/null || exit 0; sleep 1; done; kill -9 $PID 2>/dev/null ) &
  WD=$!
  wait $PID
  RC=$?
  S1=$(date +%s.%N)
  wait $WD 2>/dev/null
  wait $SP
  [ -f "$ROWS.prog" ] && { cat "$ROWS.prog" >> "$ROWS"; rm -f "$ROWS.prog"; }
  printf 'R\tE9\t%s\t%s\tprocess_wall_s\t%s\ts\tfrom the launch of dotnet to its exit, including the runtime start (the program rows give the phase alone)\n' "$L" "$R" "$(echo "$S1 - $S0" | bc -l)" >> "$ROWS"
  printf 'R\tE9\t%s\t%s\texit_code\t%s\t\tprogram exit code (137 = killed by the %s s watchdog)\n' "$L" "$R" "$RC" "$TIMEOUT" >> "$ROWS"
  if [ -n "$REF" ]; then if cmp -s "$OUT" "$REF"; then V=yes; else V=NO; fi; printf 'R\tE9\t%s\t%s\tcmp_output\t%s\t\tcmp %s %s\n' "$L" "$R" "$V" "$OUT" "$REF" >> "$ROWS"; fi
  printf 'R\tE9\t%s\t%s\tscratch_left\t%s\t\tentries of %s after the run\n' "$L" "$R" "$( [ -e "$SCR" ] && ls -A "$SCR" | wc -l || echo 0)" "$SCR" >> "$ROWS"
  printf 'R\tE9\t%s\t%s\trun_done\tyes\t\tall rows of this run are written\n' "$L" "$R" >> "$ROWS"
  echo "$L run $R rc=$RC $(tail -n 120 "$ROWS" | awk -F'\t' -v l="$L" -v r="$R" '$3==l && $4==r && ($5 ~ /_s$/ && $5 !~ /ns/ || $5=="ext_rssanon_peak" || $5=="cmp_output") {printf "%s=%.5g ", $5, $6}')"
  rm -rf "$SCR" "$OUT" $PIDF
}

RUNS=${E9_RUNS:-3}
SCRATCH=scratch.e9
case $GROUP in
psi)    for T in "$@"; do for i in $(seq 1 $RUNS); do run $T "$T.e9.psi" $i $SCRATCH $T.ref.psi out.$T.psi e9-psi --bwt $T.bwt --out out.$T.psi --pk-meta $T.pk.meta; done; done ;;
psi256) for T in "$@"; do for i in $(seq 1 $RUNS); do run $T "$T.e9.psi256" $i $SCRATCH $T.ref.psi out.$T.psi e9-psi --bwt $T.bwt --out out.$T.psi --symbols-per-pass 256; done; done ;;
bwt)    for T in "$@"; do for i in $(seq 1 $RUNS); do run $T "$T.e9.bwt" $i $SCRATCH $T.bwt out.$T.bwt e9-bwt --text $T.pk --sa $T.sa --out out.$T.bwt; done; done ;;
isa)    for T in "$@"; do
          N=$(awk -F= '$1=="n" {print $2}' $T.pk.meta); KB=$(( 64 * N / 100000000 )); [ $KB -lt 1 ] && KB=1; [ $KB -gt 64 ] && KB=64   # the chunk write buffer scaled like R: 64 KB at 10^8 (3 % of a chunk's 2.1 MB of tuples)
          for i in $(seq 1 $RUNS); do run $T "$T.e9.isa" $i $SCRATCH $T.ref.isa out.$T.isa e9-isa --sa $T.sa --out out.$T.isa --scratch $SCRATCH --scale-n 3000000000 --io-kb $KB; done
        done ;;
psiisa) for T in "$@"; do for i in $(seq 1 $RUNS); do run $T "$T.e9.psiisa" $i $SCRATCH $T.ref.psi out.$T.psi e9-psi-isa --sa $T.sa --isa $T.ref.isa --out out.$T.psi; done; done ;;
base)   if [ ! -f tiny.bwt ]; then
          dotnet "$DLL" pack --family dna:1000 --out tiny.pk > /dev/null && dotnet "$DLL" oracle --family dna:1000 --out tiny.sa > /dev/null && dotnet "$DLL" e9-bwt --text tiny.pk --sa tiny.sa --out tiny.bwt > /dev/null
        fi
        for i in $(seq 1 $RUNS); do run tiny "e9.base" $i $SCRATCH "" out.tiny.psi e9-psi --bwt tiny.bwt --out out.tiny.psi; done ;;
*) echo "unknown group $GROUP"; exit 2 ;;
esac
echo protocol-done
