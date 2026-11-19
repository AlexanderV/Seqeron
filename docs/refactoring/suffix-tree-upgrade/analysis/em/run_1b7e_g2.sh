#!/bin/bash
# Gate G2 of work package 1b.7e at 10^8 symbols (outside the shipped code): the whole container build SesaBuilder.Build (PH0 -> PH1-PH3 -> PH4-PH5 -> PH6 -> seal, one process, 4 threads, the E6 parameters
# ScaleN 3x10^9, with PSI, cold cache) on the texts of the prototype's e6 (/tmp/stbench/e6/dna_100m.txt, /tmp/stbench/e7/genome_100m.txt). Then, with decoders written outside the product (sesa_sections.py reads the header
# and the table from the text of ADR-003, crc32c_ref.c is a table CRC-32C): every section CRC-32C is recomputed, the SHA-256 of the payloads of TEXT (the .pk), SA, LCP, LCP_OVERFLOW (absent = empty) and CHILD are
# compared with the reference hashes of E7 (/tmp/stbench/e7/<name>.refs.sha), PSI with the reference of E9 (/tmp/stbench/e9/<name>.refs.sha, or, when E9 made none, derived here with the prototype of E9 from the
# E6 files: e9-bwt + e9-psi, and cross-checked with the in-RAM oracle e9-oracle), the header with the statistics of e6, and e8-verify (EmProto, eleven checks) is run on the sections extracted to files.
#   usage: run_1b7e_g2.sh <rows file> <names: dna_100m genome_100m>
set -u
ROWS=$(readlink -f "$1"); shift
NAMES=${*:-"dna_100m genome_100m"}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_1b7c_proto/SuffixTree.EmProto.dll}
E9DLL=${E9_DLL:-/tmp/stbench/bin_e9_run/SuffixTree.EmProto.dll}
D=${G2_DIR:-/tmp/stbench/c7}
W=${G2_WORK:-/tmp/stbench/g2e}
CRC=${CRC_BIN:-/tmp/stbench/crc32c_ref}
[ -x "$CRC" ] || gcc -O2 -o "$CRC" "$HERE/crc32c_ref.c"
mkdir -p "$W"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" > "$ROWS"
printf '#label\tvariant\tscan_s\tcreate_s\tpack_s\tph1_s\tph23_s\tph45_s\tbwt_s\tpsi_s\tcrc_s\tseal_s\ttotal_s\tread_bytes\twrite_bytes\tmapped_bytes\tscratch_peak\tarena_peak\trss_ph0_mb\trss_ph1_mb\trss_ph23_mb\trss_ph45_mb\trss_ph6_mb\trss_seal_mb\trss_whole_mb\tverdict\tfile_size\tgrown\treservation\toverflow\tlcpkind\trssanon_before_kb\trssanon_max_kb\trss_samples\n' >> "$ROWS"
for NAME in $NAMES; do
  rm -rf "$W/$NAME".* "$W/x".*
  case $NAME in dna_100m) TXT=/tmp/stbench/e6/dna_100m.txt ;; *) TXT=/tmp/stbench/e7/$NAME.txt ;; esac
  sync; echo 3 > /proc/sys/vm/drop_caches
  python3 -I "$HERE/sample_rss_1b7a.py" "$ROWS" "$NAME.g2" container &
  SP=$!
  G3_VARIANT=container G3_INPUT=$TXT G3_TARGET=$W/$NAME.sesa G3_ROWS=$ROWS G3_LABEL="$NAME.g2" G3_SCALE=3000000000 G3_PSI=1 \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaBuildCostTests" > "$W/last.log" 2>&1 || { echo "run failed: $NAME"; tail -20 "$W/last.log"; }
  wait $SP
  python3 -I "$HERE/sesa_sections.py" "$W/$NAME.sesa" hash --crc-c "$CRC" > "$W/$NAME.hash"
  grep -c "crc ok" "$W/$NAME.hash" | sed "s/^/#crc\t$NAME\tsections-with-crc-ok\t/" >> "$ROWS"
  grep -c "CRC MISMATCH" "$W/$NAME.hash" | sed "s/^/#crc\t$NAME\tsections-with-crc-mismatch\t/" >> "$ROWS"
  for PAIR in pk:TEXT sa:SA lcp:LCP cld:CHILD ovf:LCP_OVERFLOW; do
    X=${PAIR%%:*}; K=${PAIR##*:}
    GOT=$(awk -v k=$K '$2==k {print $5}' "$W/$NAME.hash"); [ -z "$GOT" ] && GOT=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
    WANT=$(awk -v k=$X '$1==k {print $2}' "/tmp/stbench/e7/$NAME.refs.sha")
    [ "$GOT" = "$WANT" ] && R=equal || R=DIFFERENT
    printf '#sha\t%s\t%s\t%s\t%s\n' "$NAME" "$K" "$R" "$GOT" >> "$ROWS"
  done
  WANTPSI=$(awk '$1=="psi" {print $2}' "/tmp/stbench/e9/$NAME.refs.sha" 2>/dev/null)
  SRC=e9-reference
  if [ -z "$WANTPSI" ]; then
    SRC=derived-with-the-e9-prototype
    dotnet "$E9DLL" e9-bwt --text "$D/$NAME.pk" --sa "$D/$NAME.sa" --out "$W/ref.bwt" --no-sampler > /dev/null 2>&1
    dotnet "$E9DLL" e9-psi --bwt "$W/ref.bwt" --out "$W/ref.psi" --window-kb 1024 --pk-meta "$D/$NAME.pk.meta" --no-sampler > /dev/null 2>&1
    WANTPSI=$(sha256sum "$W/ref.psi" | cut -d' ' -f1)
    rm -f "$W/ref.bwt" "$W/ref.bwt.meta" "$W/ref.psi"
  fi
  GOT=$(awk '$2=="PSI" {print $5}' "$W/$NAME.hash")
  [ "$GOT" = "$WANTPSI" ] && R=equal || R=DIFFERENT
  printf '#sha\t%s\tPSI\t%s\t%s\t%s\n' "$NAME" "$R" "$GOT" "$SRC" >> "$ROWS"
  python3 -I "$HERE/sesa_sections.py" "$W/$NAME.sesa" extract "$W/x" --no-psi
  cmp -s "$W/x.stats" "$D/$NAME.stats" && printf '#sha\t%s\tstats\tequal-to-e6-stats\t\n' "$NAME" >> "$ROWS" || printf '#sha\t%s\tstats\tDIFFERENT\t\n' "$NAME" >> "$ROWS"
  rm -f "$W/$NAME.sesa"
  dotnet "$DLL" e8-verify --prefix "$W/x" --input "$TXT" --threads 4 --scratch "$W" --case "$NAME.g2" --rows "$ROWS" > "$W/verify_$NAME.log" 2>&1
  CODE=$?
  printf '#verify\t%s\texit=%s\t%s\n' "$NAME" "$CODE" "$(grep -F "$NAME.g2" "$ROWS" | grep -F verify_all | tail -1 | cut -f5,7)" >> "$ROWS"
  rm -rf "$W/$NAME".* "$W/x".* "$W/scratch"
done
echo g2-done
