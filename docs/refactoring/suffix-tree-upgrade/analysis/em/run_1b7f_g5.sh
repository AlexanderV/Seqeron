#!/bin/bash
# Gate G5 of work package 1b.7f (outside the shipped code): the dial of the structure level, the chunk of text positions that one scan of the suffix array serves in the Phi pass (SesaVerifierTuning.ChunkSymbols, default 4 Mi).
# A container of dna_100m is built by SesaBuilder (harness variant `container`, the E6 geometry scaled to 3x10^9, PSI on, 4 threads) and opened at the Structure level (Checked) once per chunk size, cold cache, one process per run:
# wall time, RssAnon before and at the peak, the buffers and the bitmap that the pass reports.
#   usage: run_1b7f_g5.sh <rows file> [chunk sizes, default "1048576 4194304 8388608 16777216"]
set -u
ROWS=$(readlink -f "$1"); shift
CHUNKS=${*:-"1048576 4194304 8388608 16777216"}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
BUILD_PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
READ_PROJ=$ROOT/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj
W=${G5_WORK:-/tmp/stbench/g5f}
mkdir -p "$W"
TXT=/tmp/stbench/e6/dna_100m.txt
FILE=$W/dna_100m.sesa
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" > "$ROWS"
printf '#container: label variant scan create pack ph1 ph23 ph45 bwt psi crc seal total read write mapped scratch_peak arena_peak rss_ph0 rss_ph1 rss_ph23 rss_ph45 rss_ph6 rss_seal rss_whole verdict file_size grown reservation overflow lcpkind rss_before_kb rss_max_kb samples\n' >> "$ROWS"
printf '#open: label level access seconds header_s sections_s structure_s crc sectionbytes bitmapbytes bufferbytes rssanon_before_kb rssanon_max_kb rssanon_after_kb n\n' >> "$ROWS"
rm -f "$FILE"
sync; echo 3 > /proc/sys/vm/drop_caches
G3_VARIANT=container G3_INPUT=$TXT G3_TARGET=$FILE G3_ROWS=$ROWS G3_LABEL="dna_100m.g5build" G3_SCALE=3000000000 G3_PSI=1 \
  dotnet test "$BUILD_PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaBuildCostTests" > "$W/build.log" 2>&1 || { echo "build failed"; tail -20 "$W/build.log"; exit 1; }
ls -l "$FILE" | awk '{print "#file\t" $9 "\t" $5}' >> "$ROWS"
for CH in $CHUNKS; do
  sync; echo 3 > /proc/sys/vm/drop_caches
  env -u DOTNET_GCHeapHardLimit SESA_ROWS=$ROWS SESA_LABEL="dna_100m.chunk$CH" SESA_VARIANT=open SESA_LEVEL=Structure SESA_ACCESS=Checked SESA_CHUNK=$CH SESA_FILE=$FILE SESA_INPUT=$TXT \
    dotnet test "$READ_PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaReaderCostTests" > "$W/last.log" 2>&1 || { echo "reader run failed: $CH"; tail -20 "$W/last.log"; }
done
rm -f "$FILE"
