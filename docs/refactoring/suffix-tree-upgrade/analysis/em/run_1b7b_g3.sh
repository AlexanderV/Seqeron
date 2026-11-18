#!/bin/bash
# Gate G3 of work package 1b.7b (outside the shipped code): cost parity of the production MappedTextReader with the prototype PackedText, of the production mapped writer with the prototype's, and the
# measurement of fallocate against the sparse creation on the write pattern of E9. One run = one process (dotnet test of the explicit classes TextCostTests / WriterCostTests, EmProto.Tests), the writer runs with a cold
# cache (sync + drop_caches) and the device writes are read from /sys/block/<dev>/stat before the run and after the sync that follows it. DOTNET_GCHeapHardLimit = 96 MiB as in the E-series.
#   usage: run_1b7b_g3.sh <rows file> <part: reader | readertc0 | lce | writer | window> [runs, default 3]
# environment: G3_PK (default /tmp/stbench/e9/dna_100m.pk), G3_BWT (/tmp/stbench/e9/dna_100m.bwt), G3_WORK (/tmp/stbench/g3b), G3_DEV (vda)
set -u
ROWS=$(readlink -f "$1")
PART=$2
RUNS=${3:-3}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
PK=${G3_PK:-/tmp/stbench/e9/dna_100m.pk}
BWT=${G3_BWT:-/tmp/stbench/e9/dna_100m.bwt}
W=${G3_WORK:-/tmp/stbench/g3b}
DEV=${G3_DEV:-vda}
mkdir -p "$W"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" >> "$ROWS"
sect() { awk '{print $7}' /sys/block/$DEV/stat; }      # sectors written
sectr() { awk '{print $3}' /sys/block/$DEV/stat; }     # sectors read

reader() {
  printf '#cols\tlabel\twhat\tvariant\tops\tmedian_ns\tmin_ns\tmax_ns\tall_ns\tchecksum\n' >> "$ROWS"
  for K in $(seq 1 "$RUNS"); do
    [ -n "${G3_SKIP_KEYS:-}" ] || G3_PK=$PK G3_N=100000000 G3_BITS=2 G3_ROWS=$ROWS G3_KEYS=4000000 G3_LABEL=dna_100m.$PART.run$K dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~TextCostTests.RandomKeysAndSymbols" > "$W/last.log" 2>&1 || { echo "failed"; tail -20 "$W/last.log"; }
    G3_PK=$PK G3_N=100000000 G3_BITS=2 G3_ROWS=$ROWS G3_LABEL=dna_100m.$PART.run$K G3_PAIRS=random dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~TextCostTests.CommonPrefixes" > "$W/last.log" 2>&1 || { echo "failed"; tail -20 "$W/last.log"; }
    for NAME in a_1m fib_1m; do
      G3_PK=/tmp/stbench/e9/$NAME.pk G3_N=1000000 G3_BITS=1 G3_ROWS=$ROWS G3_LABEL=$NAME.$PART.run$K G3_PAIRS=near G3_KEYS=1000000 dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~TextCostTests.CommonPrefixes" > "$W/last.log" 2>&1 || { echo "failed $NAME"; tail -20 "$W/last.log"; }
    done
  done
}

run_writer() { # label variant pattern window run
  L=$1; V=$2; P=$3; WIN=$4; K=$5
  rm -rf "$W/out"; sync
  echo 3 > /proc/sys/vm/drop_caches
  cat "$BWT" > /dev/null   # the input stays warm: only the output is cold
  S0=$(sect); R0=$(sectr); T0=$(date +%s.%N)
  G3W_VARIANT=$V G3W_PATTERN=$P G3W_BWT=$BWT G3W_OUT=$W/out G3W_WINDOW=$WIN G3W_ROWS=$ROWS G3W_LABEL="$L.run$K" \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~WriterCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $L $V $P"; tail -20 "$W/last.log"; }
  T1=$(date +%s.%N); sync; T2=$(date +%s.%N)
  S1=$(sect); R1=$(sectr)
  printf '#dev\t%s.run%s\t%s\t%s\t%s\tdevice_written_bytes\t%s\tdevice_read_bytes\t%s\twall_s\t%s\twall_with_sync_s\t%s\n' "$L" "$K" "$V" "$P" "$WIN" "$(( (S1 - S0) * 512 ))" "$(( (R1 - R0) * 512 ))" "$(echo "$T1 - $T0" | bc)" "$(echo "$T2 - $T0" | bc)" >> "$ROWS"
  du -sk "$W/out" 2>/dev/null | awk -v l="$L.run$K" '{print "#du\t" l "\tallocated_kb\t" $1}' >> "$ROWS"
  rm -rf "$W/out"
}

writer() {
  printf '#cols\tlabel\tvariant\tpattern\twindow\tseconds_create\tseconds_total\trssanon_before_kb\trssanon_max_kb\tbytes\twindows\n' >> "$ROWS"
  for P in seq psi; do
    for K in $(seq 1 "$RUNS"); do
      case $((K % 3)) in 1) ORDER="proto sparse fallocate" ;; 2) ORDER="sparse fallocate proto" ;; 0) ORDER="fallocate proto sparse" ;; esac
      for V in $ORDER; do run_writer "w" "$V" "$P" 1048576 "$K"; done
    done
  done
}

window() { # the harness check against E9: 64 KB windows of the Psi pattern
  printf '#cols\tlabel\tvariant\tpattern\twindow\tseconds_create\tseconds_total\trssanon_before_kb\trssanon_max_kb\tbytes\twindows\n' >> "$ROWS"
  for K in 1 2; do for V in proto sparse; do run_writer "w64k" "$V" psi 65536 "$K"; done; done
}

case $PART in lce) G3_SKIP_KEYS=1 reader ;; reader) reader ;; readertc0) export DOTNET_TieredCompilation=0; reader ;; writer) writer ;; window) window ;; *) echo "part: reader | readertc0 | writer | window"; exit 2 ;; esac
echo "g3-$PART-done"
