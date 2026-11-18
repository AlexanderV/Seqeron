#!/bin/bash
# Gate G2 of work package 1b.7c at 10^8 symbols (outside the shipped code): the production stage on the packed text and the suffix array of the prototype's e6 (/tmp/stbench/c7: dna_100m and genome_100m, the E7 setup, the
# suffix arrays and packed texts equal to the E7 reference hashes), written to files of its own (hard-linked inputs, the production .lcp .ovf .cld .stats): SHA-256 of every output against the reference hashes of E7
# (/tmp/stbench/e7/<name>.refs.sha) and the verifier of E8 (e8-verify: eleven checks, with the input text) on the production outputs. One process per input, 4 threads, the parameters of E6 (scale-n 3x10^9).
#   usage: run_1b7c_g2.sh <rows file> <names: dna_100m genome_100m>
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
for NAME in $NAMES; do
  rm -rf "$W/$NAME".* "$W/scratch"
  for X in pk pk.meta sa; do ln "$D/$NAME.$X" "$W/$NAME.$X"; done
  N=$(grep '^n=' "$D/$NAME.pk.meta" | cut -d= -f2); BITS=$(grep '^bits=' "$D/$NAME.pk.meta" | cut -d= -f2); SIGMA=$(grep '^sigma=' "$D/$NAME.pk.meta" | cut -d= -f2)
  case $NAME in dna_100m) TXT=/tmp/stbench/e6/dna_100m.txt ;; *) TXT=/tmp/stbench/e7/$NAME.txt ;; esac
  sync; echo 3 > /proc/sys/vm/drop_caches
  G3_VARIANT=prod G3_PK=$W/$NAME.pk G3_N=$N G3_BITS=$BITS G3_SA=$W/$NAME.sa G3_OUT=$W/$NAME G3_SCRATCH=$W/scratch G3_SCALE=3000000000 G3_SIGMA=$SIGMA G3_ROWS=$ROWS G3_LABEL="$NAME.g2" \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~LcpChildCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $NAME"; tail -20 "$W/last.log"; }
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
