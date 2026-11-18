#!/bin/bash
# Gate G3 of work package 1b.7a (outside the shipped code): cost parity of the production Dist with the prototype IsaByDist on ISA of dna_100m.
# One run = one process (dotnet test of the explicit class DistCostTests, which times the Build call alone and samples RssAnon every 10 ms with the production sampler), cold cache (sync + drop_caches) or warm cache
# (the suffix-array file read once before the run), DOTNET_GCHeapHardLimit = 96 MiB as in the E-series, variants alternate within a run number so that neither always goes first.
#   usage: run_1b7a_g3.sh <rows file> [runs per cell, default 3]
# environment: G3_SA (default /tmp/stbench/e9/dna_100m.sa), G3_N (100000000), G3_REF (the in-RAM ISA of the E9 setup: every output is cmp'ed with it outside the timed region), G3_WORK (/tmp/stbench/g3)
set -u
ROWS=$(readlink -f "$1")
RUNS=${2:-3}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
SA=${G3_SA:-/tmp/stbench/e9/dna_100m.sa}
N=${G3_N:-100000000}
REF=${G3_REF:-/tmp/stbench/e9/dna_100m.ref.isa}
W=${G3_WORK:-/tmp/stbench/g3}
mkdir -p "$W"
printf '#label\tvariant\tR\tio_buffer\tseconds\trssanon_before_kb\trssanon_max_kb\trss_samples\tarena_peak_bytes\tread_bytes\twrite_bytes\tmapped_bytes\tscratch_peak_bytes\n' > "$ROWS"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" >> "$ROWS"

run() { # label variant R io cache run
  L=$1; V=$2; R=$3; IO=$4; CACHE=$5; K=$6
  rm -rf "$W/scratch" "$W/isa.out"
  sync   # the dirty pages of the previous run (400 MB of output, 800 MB of scratch) must not be written back during this one
  if [ "$CACHE" = cold ]; then echo 3 > /proc/sys/vm/drop_caches; else cat "$SA" > /dev/null; fi
  python3 -I "$HERE/sample_rss_1b7a.py" "$ROWS" "$L.$CACHE.$K" "$V" &
  SP=$!
  G3_VARIANT=$V G3_SA=$SA G3_N=$N G3_R=$R G3_IO=$IO G3_OUT=$W/isa.out G3_SCRATCH=$W/scratch G3_ROWS=$ROWS G3_LABEL="$L.$CACHE.$K" \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~DistCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $L $V $CACHE $K"; tail -20 "$W/last.log"; }
  wait $SP
  if [ -f "$REF" ]; then if cmp -s "$W/isa.out" "$REF"; then echo "#cmp	$L.$CACHE.$K	$V	equal-to-reference" >> "$ROWS"; else echo "#cmp	$L.$CACHE.$K	$V	DIFFERENT" >> "$ROWS"; fi; fi
  rm -rf "$W/scratch" "$W/isa.out"
}

# A: the 3x10^9 geometry scaled to 10^8 (R = 266,667, P = 375, 64 KB buffers) as in E9; B: the real R = 8,000,000 (P = 13)
for CELL in "A 266667 65536" "B 8000000 65536"; do
  set -- $CELL
  for CACHE in cold warm; do
    for K in $(seq 1 "$RUNS"); do
      if [ $((K % 2)) = 1 ]; then ORDER="proto prod"; else ORDER="prod proto"; fi
      for V in $ORDER; do run "$1" "$V" "$2" "$3" "$CACHE" "$K"; done
    done
  done
done
# the counter overhead: production with and without counters, warm, parameters A, alternating
for K in $(seq 1 "$RUNS"); do
  if [ $((K % 2)) = 1 ]; then ORDER="prod prodnc"; else ORDER="prodnc prod"; fi
  for V in $ORDER; do run "Aoverhead" "$V" 266667 65536 warm "$K"; done
done
echo g3-done
