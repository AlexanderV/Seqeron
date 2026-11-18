#!/bin/bash
# Gate G2 of work package 1b.7d at 10^8 symbols (outside the shipped code): the whole production chain PH0 (PackedTextWriter) -> PH1 + PH2/PH3 (SaBuilder) -> PH4 + PH5 (LcpChildBuilder), one process, one shared IoCounters, on
# the texts of the prototype's e6 (/tmp/stbench/e6/dna_100m.txt, /tmp/stbench/e7/genome_100m.txt), with the parameters of E6 (ScaleN 3x10^9): SHA-256 of .pk .sa .lcp .ovf .cld against the reference hashes of E7
# (/tmp/stbench/e7/<name>.refs.sha), the statistics against the e6 statistics, and the verifier of E8 (e8-verify, eleven checks, with the input text) on the production outputs. 4 threads, cold cache.
#   usage: run_1b7d_g2.sh <rows file> <names: dna_100m genome_100m>
set -u
ROWS=$(readlink -f "$1"); shift
NAMES=${*:-"dna_100m genome_100m"}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_1b7c_proto/SuffixTree.EmProto.dll}
D=${G2_DIR:-/tmp/stbench/c7}
W=${G2_WORK:-/tmp/stbench/g2c}
mkdir -p "$W"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" > "$ROWS"
printf '#label\tvariant\tph0_s\tph1_s\tph23_s\tph45_s\trssanon_before_kb\trssanon_max_kb\trss_samples\tsa_rss_mb\tsa_verdict\tlcp_verdict\tread_bytes\twrite_bytes\tmapped_bytes\tscratch_peak\n' >> "$ROWS"
for NAME in $NAMES; do
  rm -rf "$W/$NAME".* "$W/scratch"
  SIGMA=$(grep '^sigma=' "$D/$NAME.pk.meta" | cut -d= -f2)
  case $NAME in dna_100m) TXT=/tmp/stbench/e6/dna_100m.txt ;; *) TXT=/tmp/stbench/e7/$NAME.txt ;; esac
  sync; echo 3 > /proc/sys/vm/drop_caches
  G3_VARIANT=e2e G2_TXT=$TXT G2_OUT=$W/$NAME G2_SIGMA=$SIGMA G3_SCRATCH=$W/scratch G3_SCALE=3000000000 G3_ROWS=$ROWS G3_LABEL="$NAME.g2" \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SaBuildCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $NAME"; tail -20 "$W/last.log"; }
  for X in pk sa lcp ovf cld; do
    GOT=$(sha256sum "$W/$NAME.$X" | cut -d' ' -f1); WANT=$(awk -v k=$X '$1==k {print $2}' "/tmp/stbench/e7/$NAME.refs.sha")
    [ "$GOT" = "$WANT" ] && R=equal || R=DIFFERENT
    printf '#sha\t%s\t%s\t%s\t%s\n' "$NAME" "$X" "$R" "$GOT" >> "$ROWS"
  done
  cmp -s "$W/$NAME.stats" "$D/$NAME.stats" && printf '#sha\t%s\tstats\tequal-to-e6-stats\t\n' "$NAME" >> "$ROWS" || printf '#sha\t%s\tstats\tDIFFERENT\t\n' "$NAME" >> "$ROWS"
  dotnet "$DLL" e8-verify --prefix "$W/$NAME" --input "$TXT" --threads 4 --scratch "$W" --case "$NAME.g2" --rows "$ROWS" > "$W/verify_$NAME.log" 2>&1
  CODE=$?
  printf '#verify\t%s\texit=%s\t%s\n' "$NAME" "$CODE" "$(grep -F "$NAME.g2" "$ROWS" | grep -F verify_all | tail -1 | cut -f5,7)" >> "$ROWS"
  rm -rf "$W/$NAME".* "$W/scratch"
done
echo g2-done
