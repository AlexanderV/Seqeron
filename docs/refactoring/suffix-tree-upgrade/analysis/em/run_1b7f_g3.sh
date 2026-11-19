#!/bin/bash
# Gate G3 / G5 of work package 1b.7f (outside the shipped code): a container of dna_10m (10 Mbp) built by SesaBuilder (harness variant `container`, real geometry), then
#   cold-open   10 runs of the Header-level open, cache dropped (sync + echo 3 > /proc/sys/vm/drop_caches) before each (one process per run: the time of the open call, mapping and header validation, not of the process start)
#   levels      Header, Sections and Structure at 10 Mbp, cold, 3 runs each
#   throughput  count of 20-mers, locate of 8-mers, count of absent 24-mers and a child walk, in memory against the file at Checked and Trusted, warm, median of 5 per process, one tree per process, 3 rounds in alternating orders
#   usage: run_1b7f_g3.sh <rows file> [runs, default 10]
set -u
ROWS=$(readlink -f "$1")
RUNS=${2:-10}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
BUILD_PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
READ_PROJ=$ROOT/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj
W=${G3_WORK:-/tmp/stbench/g3f}
mkdir -p "$W"
TXT=/tmp/stbench/data/dna_10m.txt
FILE=$W/dna_10m.sesa
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" > "$ROWS"
printf '#container: label variant scan create pack ph1 ph23 ph45 bwt psi crc seal total read write mapped scratch_peak arena_peak rss_ph0 rss_ph1 rss_ph23 rss_ph45 rss_ph6 rss_seal rss_whole verdict file_size grown reservation overflow lcpkind rss_before_kb rss_max_kb samples\n' >> "$ROWS"
printf '#open: label level access seconds header_s sections_s structure_s crc sectionbytes bitmapbytes bufferbytes rssanon_before_kb rssanon_max_kb rssanon_after_kb n\n' >> "$ROWS"
printf '#throughput: label tree count20_per_s locate8_per_s count24absent_per_s childwalk_nodes_per_s repeats\n' >> "$ROWS"
rm -f "$FILE"
G3_VARIANT=container G3_INPUT=$TXT G3_TARGET=$FILE G3_ROWS=$ROWS G3_LABEL="dna_10m.build" G3_SCALE=0 G3_PSI=1 \
  dotnet test "$BUILD_PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaBuildCostTests" > "$W/build.log" 2>&1 || { echo "build failed"; tail -20 "$W/build.log"; exit 1; }
ls -l "$FILE" | awk '{print "#file\t" $9 "\t" $5}' >> "$ROWS"
reader() { # label level access [extra env]
  sync; echo 3 > /proc/sys/vm/drop_caches
  env -u DOTNET_GCHeapHardLimit SESA_ROWS=$ROWS ${EXTRA:-} SESA_LABEL="$1" SESA_VARIANT=open SESA_LEVEL=$2 SESA_ACCESS=${3:-Checked} SESA_FILE=$FILE SESA_INPUT=$TXT \
    dotnet test "$READ_PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaReaderCostTests" > "$W/last.log" 2>&1 || { echo "reader run failed: $1 $2"; tail -20 "$W/last.log"; }
}
for K in $(seq 1 "$RUNS"); do EXTRA="SESA_WARM=1" reader "cold-open.$K" Header Checked; done
for K in 1 2 3; do reader "levels.$K" Sections Checked; reader "levels.$K" Structure Checked; done
# one tree per process, three rounds in alternating orders (a process that holds the in-memory tree first and the others after it would carry the order as a bias)
ROUND=0
for ORDER in "memory checked trusted" "trusted checked memory" "checked memory trusted"; do
  ROUND=$((ROUND + 1))
  for TREE in $ORDER; do
    sync
    env -u DOTNET_GCHeapHardLimit SESA_ROWS=$ROWS SESA_LABEL=dna_10m.r$ROUND SESA_VARIANT=throughput SESA_TREES=$TREE SESA_FILE=$FILE SESA_INPUT=$TXT SESA_REPEATS=5 \
      dotnet test "$READ_PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaReaderCostTests" > "$W/throughput.log" 2>&1 || { echo "throughput failed: $TREE"; tail -20 "$W/throughput.log"; }
  done
done
