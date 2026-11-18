#!/bin/bash
# Gate G3 of work package 1b.7c (outside the shipped code): cost parity of the production PH4 + PH5 (LcpChildBuilder) with the prototype's fused PH4 + PH5 (ExternalLcp + E5.LazyScan, the parameters of ScaleRules) on a given
# packed text and suffix array file. One run = one process (dotnet test of the explicit class LcpChildCostTests, which times the build call alone and samples RssAnon every 10 ms with the production sampler; an independent
# sampler, sample_rss_1b7a.py, polls the test host from outside), cold cache (sync + drop_caches) or warm cache (the SA and the packed text read once before the run), DOTNET_GCHeapHardLimit = 96 MiB as in the E-series, the
# variants alternate within a run number so that neither always goes first. After every run the outputs are cmp'ed with the prototype's reference outputs outside the timed region and deleted.
#   usage: run_1b7c_g3.sh <rows file> <inputs: dna genome an> [runs per cell, default 3]
# environment: G3_DIR (default /tmp/stbench/c7: <name>.pk, <name>.sa and the reference outputs <name>.lcp .ovf .cld of the prototype's e6), G3_WORK (/tmp/stbench/g3c)
set -u
ROWS=$(readlink -f "$1")
INPUTS=${2:-"dna genome"}
RUNS=${3:-3}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
D=${G3_DIR:-/tmp/stbench/c7}
W=${G3_WORK:-/tmp/stbench/g3c}
mkdir -p "$W"
printf '#label\tvariant\tseconds\trssanon_before_kb\trssanon_max_kb\trss_samples\tarena_peak_bytes\tchild_arena_peak_bytes\tread_bytes\twrite_bytes\tmapped_lcp_bytes\tchild_write_bytes\tspill_written\tspill_read\tchild_flushed\tchild_reopened\tstack_peak\tb2_workers\tR\tscratch_peak\tb1_s\tb2_s\tb3_s\tcompare_steps\n' > "$ROWS"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" >> "$ROWS"

meta() { grep "^$2=" "$1.pk.meta" | cut -d= -f2; }

run() { # name label variant scale cache run
  NAME=$1; L=$2; V=$3; SC=$4; CACHE=$5; K=$6
  N=$(meta "$D/$NAME" n); BITS=$(meta "$D/$NAME" bits)
  rm -rf "$W/scratch" "$W/out".*
  sync   # the dirty pages of the previous run (500 MB of output) must not be written back during this one
  if [ "$CACHE" = cold ]; then echo 3 > /proc/sys/vm/drop_caches; else cat "$D/$NAME.sa" "$D/$NAME.pk" > /dev/null; fi
  python3 -I "$HERE/sample_rss_1b7a.py" "$ROWS" "$L.$CACHE.$K" "$V" &
  SP=$!
  G3_VARIANT=$V G3_PK=$D/$NAME.pk G3_N=$N G3_BITS=$BITS G3_SA=$D/$NAME.sa G3_OUT=$W/out G3_SCRATCH=$W/scratch G3_SCALE=$SC G3_ROWS=$ROWS G3_LABEL="$L.$CACHE.$K" \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~LcpChildCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $L $V $CACHE $K"; tail -20 "$W/last.log"; }
  wait $SP
  EQ=equal-to-reference
  for X in lcp ovf cld; do cmp -s "$W/out.$X" "$D/$NAME.$X" || EQ="DIFFERENT($X)"; done
  echo "#cmp	$L.$CACHE.$K	$V	$EQ" >> "$ROWS"
  rm -rf "$W/scratch" "$W/out".*
}

for NAME in $INPUTS; do
  # A: the 3x10^9 geometry scaled to 10^8 (R = 266,667, P = 375) as in E6 / E7; B: the real R = 8,000,000 (P = 13); the cold cell of A is the gate, the others are reported
  for CELL in "A 3000000000 cold" "A 3000000000 warm" "B 0 cold"; do
    set -- $CELL
    for K in $(seq 1 "$RUNS"); do
      if [ $((K % 2)) = 1 ]; then ORDER="proto prod"; else ORDER="prod proto"; fi
      for V in $ORDER; do run "$NAME" "$NAME.$1" "$V" "$2" "$3" "$K"; done
    done
  done
  # the cost of the builder's own RSS sampler: production with and without it, warm, parameters A
  for K in $(seq 1 "$RUNS"); do
    if [ $((K % 2)) = 1 ]; then ORDER="prod prodnr"; else ORDER="prodnr prod"; fi
    for V in $ORDER; do run "$NAME" "$NAME.Aoverhead" "$V" 3000000000 warm "$K"; done
  done
done
echo g3-done
