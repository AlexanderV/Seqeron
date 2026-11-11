#!/bin/bash
# E4 real-size runs (outside the shipped code): UNSCALED parameters (R = 8,000,000, tuple buffers 32 KiB, 160 MiB arena cap, B2 tasks from the reference geometry), cold cache, 4 threads,
# DOTNET_GCHeapHardLimit = 96 MiB. The text is packed, the suffix array comes from the E3 builder (its phase is timed), the oracle is the in-RAM LcpKernel where it fits (random DNA) or the structure and
# sampled-rank check of `check-lcp` (a^n: the kernel would need about 6.5 GB for n = 500 M).
# usage: run_e4_big.sh <work dir> <rows file> <name> <family: rand:N | a:N | genome:N>
set -u
W=$1; ROWS=$2; X=$3; FAM=$4
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e4_run/SuffixTree.EmProto.dll}
mkdir -p "$W/scratch"; cd "$W" || exit 1
printf 'R\tE4\tbig\thead\t%s\t\tgit rev-parse HEAD of the binary under test; %s %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$X" "$FAM" >> "$ROWS"
export DOTNET_GCHeapHardLimit=0x6000000
if [ "${SKIP_SETUP:-no}" != yes ]; then
  case $FAM in
    rand:*) N=${FAM#rand:}; ( unset DOTNET_GCHeapHardLimit; dotnet "$DLL" packrand --n "$N" --out "$X.pk" > /dev/null ) ;;
    *) ( unset DOTNET_GCHeapHardLimit; dotnet "$DLL" pack --family "$FAM" --out "$X.pk" > /dev/null ) ;;
  esac
  ( unset DOTNET_GCHeapHardLimit; dotnet "$DLL" mkranks --text "$X.pk" --out "$X.ranks" > /dev/null )
  sync; echo 3 > /proc/sys/vm/drop_caches
  dotnet "$DLL" e3 --text "$X.pk" --ranks "$X.ranks" --sa "$X.sa" --scratch scratch --threads 4 --case "$X.sa_phase" --run 1 --rows "$ROWS" > /dev/null
  echo "$X SA rc=$?"
fi
HAVE_ORACLE=no
if [ "${ORACLE:-yes}" = yes ]; then
  ( unset DOTNET_GCHeapHardLimit; dotnet "$DLL" lcp-oracle --text "$X.pk" --sa "$X.sa" --out "$X.oracle" --case "$X.kernel" --rows "$ROWS" > /dev/null ) && HAVE_ORACLE=yes
  echo "$X oracle rc=$? have=$HAVE_ORACLE"
fi
run() { # label run extra-args...
  L=$1; R=$2; shift 2
  V=""; [ "$HAVE_ORACLE" = yes ] && V="--verify $X.oracle"
  sync; echo 3 > /proc/sys/vm/drop_caches
  dotnet "$DLL" e4 --text "$X.pk" --sa "$X.sa" --out "$X.out" --scratch scratch $V --case "$L" --run "$R" --rows "$ROWS" "$@" > /dev/null
  echo "$L run $R rc=$?"
  if [ -n "${ALLSAME_N:-}" ]; then # exact closed form of a^n for every run (no oracle fits in RAM)
    if python3 -I "$(dirname "$(readlink -f "$0")")/check_allsame.py" "$ALLSAME_N" "$X.sa" "$X.out.lcp" "$X.out.ovf" > /dev/null; then Y=yes; else Y=NO; fi
    printf 'R\tE4\t%s\t%s\tlcp_equals_oracle\t%s\t\tbase and overflow equal the closed form of a^n for every rank (check_allsame.py)\n' "$L" "$R" "$Y" >> "$ROWS"
  fi
}
for i in 1 2 3; do run "$X.dist.t4" $i --threads 4; done
if [ "$HAVE_ORACLE" != yes ] && [ -z "${ALLSAME_N:-}" ]; then
  dotnet "$DLL" check-lcp --text "$X.pk" --sa "$X.sa" --lcp "$X.out" --samples 2000 --case "$X.dist.check" --run 1 --rows "$ROWS" > /dev/null; echo "$X check rc=$?"
  cp "$X.out.lcp" "$X.dist.lcp.keep"; cp "$X.out.ovf" "$X.dist.ovf.keep"
fi
for i in 1 2 3; do run "$X.b2.t4" $i --threads 4 --mode bounded --rounds 2; done
if [ "$HAVE_ORACLE" != yes ] && [ -z "${ALLSAME_N:-}" ]; then
  cmp "$X.out.lcp" "$X.dist.lcp.keep" && cmp "$X.out.ovf" "$X.dist.ovf.keep" && printf 'R\tE4\t%s\t1\tbounded_equals_distributed\tyes\t\tbyte-identical base and overflow files of the two modes (cmp)\n' "$X.b2.t4" >> "$ROWS"
  rm -f "$X.dist.lcp.keep" "$X.dist.ovf.keep"
fi
run "$X.b1.t4" 1 --threads 4 --mode bounded --rounds 1
run "$X.b4.t4" 1 --threads 4 --mode bounded --rounds 4
run "$X.dist.t1" 1 --threads 1
run "$X.dist.t2" 1 --threads 2
run "$X.b2.t1" 1 --threads 1 --mode bounded --rounds 2
# the arena cap enforced from inside: four B2 tasks do not fit 160 MiB at R = 8,000,000
sync; echo 3 > /proc/sys/vm/drop_caches
dotnet "$DLL" e4 --text "$X.pk" --sa "$X.sa" --out "$X.out" --scratch scratch --threads 4 --b2-workers 4 --case "$X.enforce" --run 1 --rows "$ROWS" > "$X.enforce.out" 2> "$X.enforce.err"
echo "$X enforce rc=$?"; head -c 600 "$X.enforce.err" > "$ROWS.enforce.txt"
echo big-done
