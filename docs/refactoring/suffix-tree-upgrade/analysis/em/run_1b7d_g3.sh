#!/bin/bash
# Gate G3 / G5 of work package 1b.7d (outside the shipped code): cost parity of the production PH1 and PH2/PH3 (SaBuilder) with the prototype's (ExternalDcs + EmBuilder, the parameters of ScaleRules, phase by phase as
# E6 runs them) on a given packed text. One run = one process (dotnet test of the explicit class SaBuildCostTests, which times the phases and samples RssAnon every 10 ms with the production sampler; an independent sampler,
# sample_rss_1b7a.py, polls the test host from outside), cold cache (sync + drop_caches) or warm cache (the packed text read once before the run), DOTNET_GCHeapHardLimit = 96 MiB as in the E-series, the variants alternate
# within a run number so that neither always goes first. After every run the suffix array is cmp'ed with the reference outside the timed region and deleted.
#   usage: run_1b7d_g3.sh <rows file> <inputs: dna genome> [runs per cell, default 3] [cells, default "A-cold B-cold A-warm"]
# environment: G3_DIR (default /tmp/stbench/c7: <name>_100m.pk and the reference <name>_100m.sa of the prototype's e6), G3_WORK (/tmp/stbench/g3c)
#   cells: A = the 3x10^9 geometry scaled to 10^8 (ScaleN 3x10^9, the E6 parameters), B = the real geometry (ScaleN 0); G5 = the cost of the reservation (prod against prodsparse) and of the sampler (prod against prodnr), A warm; G5v = the cost of the optional rank-file check (prod against prodverify), A warm
set -u
ROWS=$(readlink -f "$1")
INPUTS=${2:-"dna genome"}
RUNS=${3:-3}
CELLS=${4:-"A-cold B-cold A-warm"}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
D=${G3_DIR:-/tmp/stbench/c7}
W=${G3_WORK:-/tmp/stbench/g3c}
mkdir -p "$W"
printf '#label\tvariant\tph1_s\tph23_s\ttotal_s\tph1_read\tph1_write\tph1_mapped\tph23_read\tph23_write\tph23_mapped\tscratch_peak_ph1\tscratch_peak\tarena_ph1\tarena_ph23\tlevels\tverdict\trssanon_before_kb\trssanon_max_kb\trss_samples\n' > "$ROWS"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" >> "$ROWS"

meta() { grep "^$2=" "$1.pk.meta" | cut -d= -f2; }

run() { # name label variant scale cache run
  NAME=$1; L=$2; V=$3; SC=$4; CACHE=$5; K=$6
  N=$(meta "$D/${NAME}_100m" n); BITS=$(meta "$D/${NAME}_100m" bits)
  rm -rf "$W/scratch" "$W/out".*
  sync   # the dirty pages of the previous run (400 MB of output) must not be written back during this one
  if [ "$CACHE" = cold ]; then echo 3 > /proc/sys/vm/drop_caches; else cat "$D/${NAME}_100m.pk" > /dev/null; fi
  python3 -I "$HERE/sample_rss_1b7a.py" "$ROWS" "$L.$CACHE.$K" "$V" &
  SP=$!
  G3_VARIANT=$V G3_PK=$D/${NAME}_100m.pk G3_N=$N G3_BITS=$BITS G3_SA=$W/out.sa G3_SCRATCH=$W/scratch G3_SCALE=$SC G3_ROWS=$ROWS G3_LABEL="$L.$CACHE.$K" \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SaBuildCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $L $V $CACHE $K"; tail -20 "$W/last.log"; }
  wait $SP
  EQ=equal-to-reference
  cmp -s "$W/out.sa" "$D/${NAME}_100m.sa" || EQ="DIFFERENT(sa)"
  echo "#cmp	$L.$CACHE.$K	$V	$EQ" >> "$ROWS"
  rm -rf "$W/scratch" "$W/out".*
}

for NAME in $INPUTS; do
  for CELL in $CELLS; do
    case $CELL in
      A-cold) SC=3000000000; CACHE=cold; VARS="proto prod" ;;
      B-cold) SC=0; CACHE=cold; VARS="proto prod" ;;
      A-warm) SC=3000000000; CACHE=warm; VARS="proto prod" ;;
      G5-A-warm) SC=3000000000; CACHE=warm; VARS="prod prodsparse prodnr" ;;
      G5v-A-warm) SC=3000000000; CACHE=warm; VARS="prod prodverify" ;;
      *) echo "unknown cell $CELL"; continue ;;
    esac
    for K in $(seq 1 "$RUNS"); do
      set -- $VARS
      # rotate the order within the cell so that no variant always goes first
      ORDER=$(for i in $(seq 0 $(($#-1))); do eval echo \${$(( (i + K - 1) % $# + 1 ))}; done)
      for V in $ORDER; do run "$NAME" "$NAME.$CELL" "$V" "$SC" "$CACHE" "$K"; done
    done
  done
done
echo g3-done
