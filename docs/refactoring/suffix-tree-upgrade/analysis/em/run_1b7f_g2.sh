#!/bin/bash
# Gate G2 of work package 1b.7f at 10^8 symbols (outside the shipped code). The container is built by SesaBuilder.Build through the harness variant `container` of 1b.7e (SesaBuildCostTests, PSI on, 4 threads, the E6 parameters
# ScaleN 3x10^9) from /tmp/stbench/e6/dna_100m.txt and /tmp/stbench/e7/genome_100m.txt, and then read by the reader of 1b.7f (SesaReaderCostTests, one variant per process, cold cache before every process):
#   open Header, open Sections, open Structure (wall time, RssAnon before and at the peak, the validation report),
#   identity (1,000 patterns of the text and 1,000 random ones: counts and locate sets against an independent binary search over the suffix array and the input text),
#   e8-verify --container through the adapter of EmProto (eleven checks, the sections read in place),
#   corrupt (an adjacent swap of the suffix array and a child cell, the section CRC recomputed: the structure level must catch both; the container is restored afterwards).
#   usage: run_1b7f_g2.sh <rows file> <names: dna_100m genome_100m>
set -u
ROWS=$(readlink -f "$1"); shift
NAMES=${*:-"dna_100m genome_100m"}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
BUILD_PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
READ_PROJ=$ROOT/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_1b7f/SuffixTree.EmProto.dll}
W=${G2_WORK:-/tmp/stbench/g2f}
mkdir -p "$W"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" > "$ROWS"
printf '#container: label variant scan create pack ph1 ph23 ph45 bwt psi crc seal total read write mapped scratch_peak arena_peak rss_ph0 rss_ph1 rss_ph23 rss_ph45 rss_ph6 rss_seal rss_whole verdict file_size grown reservation overflow lcpkind rss_before_kb rss_max_kb samples\n' >> "$ROWS"
printf '#open: label level access seconds header_s sections_s structure_s crc sectionbytes bitmapbytes bufferbytes rssanon_before_kb rssanon_max_kb rssanon_after_kb n\n' >> "$ROWS"
printf '#identity: label level patterns absent located_fully located_prefix count_mismatches locate_mismatches seconds\n' >> "$ROWS"
printf '#corrupt: label case verdict seconds message\n' >> "$ROWS"

reader() { # label variant level access [extra env]
  sync; echo 3 > /proc/sys/vm/drop_caches
  env -u DOTNET_GCHeapHardLimit SESA_ROWS=$ROWS SESA_LABEL="$1" SESA_VARIANT=$2 SESA_LEVEL=$3 SESA_ACCESS=${4:-Checked} SESA_FILE=$FILE SESA_INPUT=$TXT ${EXTRA:-} \
    dotnet test "$READ_PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaReaderCostTests" > "$W/last.log" 2>&1 || { echo "reader run failed: $1 $2 $3"; tail -20 "$W/last.log"; }
}

for NAME in $NAMES; do
  rm -f "$W/$NAME.sesa"
  case $NAME in dna_100m) TXT=/tmp/stbench/e6/dna_100m.txt ;; *) TXT=/tmp/stbench/e7/$NAME.txt ;; esac
  FILE=$W/$NAME.sesa
  sync; echo 3 > /proc/sys/vm/drop_caches
  G3_VARIANT=container G3_INPUT=$TXT G3_TARGET=$FILE G3_ROWS=$ROWS G3_LABEL="$NAME.g2build" G3_SCALE=3000000000 G3_PSI=1 \
    dotnet test "$BUILD_PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaBuildCostTests" > "$W/build.log" 2>&1 || { echo "build failed: $NAME"; tail -20 "$W/build.log"; continue; }
  ls -l "$FILE" | awk '{print "#file\t" $9 "\t" $5}' >> "$ROWS"
  reader "$NAME.g2" open Header Checked
  reader "$NAME.g2" open Sections Checked
  reader "$NAME.g2" open Structure Checked
  EXTRA="SESA_SHA=1" reader "$NAME.g2sha" open Structure Checked
  reader "$NAME.g2" identity Header Checked
  sync; echo 3 > /proc/sys/vm/drop_caches
  dotnet "$DLL" e8-verify --container "$FILE" --input "$TXT" --threads 4 --scratch "$W" --case "$NAME.g2" --rows "$ROWS" > "$W/verify_$NAME.log" 2>&1
  CODE=$?
  printf '#e8-verify\t%s\texit\t%s\n' "$NAME" "$CODE" >> "$ROWS"
  tail -3 "$W/verify_$NAME.log" >> "$W/verify_tail_$NAME.log"
  reader "$NAME.g2" corrupt Structure Checked
  rm -f "$FILE"
done
