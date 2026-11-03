#!/usr/bin/env bash
# BWA (index -a bwtsw) and Bowtie2 (bowtie2-build): FM-index build time + index size only (no count/locate API comparable here).
# Needs: either the distro packages `bwa` `bowtie2` on PATH, or git+make+zlib to build from source.
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BWA_URL=https://github.com/lh3/bwa;                BWA_SHA=b92993c1161e73167181558856567ef2f367e3f0   # v0.7.19
BT2_URL=https://github.com/BenLangmead/bowtie2;    BT2_SHA=0c6a1c75e047ad8bf70c178fa3cb1528fba6adc2   # v2.5.5
BWA="${BWA_BIN:-}"; BT2B="${BT2_BUILD_BIN:-}"
# no binaries given -> build both from the pinned upstream commits (needs git make g++ zlib1g-dev perl)
if [ -z "$BWA" ]; then
  need git; need make; need gcc
  clone_at $BWA_URL $BWA_SHA "$WORK/bwa"; note_licence "$WORK/bwa"
  make -C "$WORK/bwa" -j2 >"$WORK/bwa_make.log" 2>&1 || { tail -20 "$WORK/bwa_make.log"; echo "bwa build failed"; }
  [ -x "$WORK/bwa/bwa" ] && BWA="$WORK/bwa/bwa"
fi
if [ -z "$BT2B" ]; then
  need git; need make; need g++
  clone_at $BT2_URL $BT2_SHA "$WORK/bowtie2"; note_licence "$WORK/bowtie2"
  git -C "$WORK/bowtie2" submodule update --init --recursive --quiet || true
  make -C "$WORK/bowtie2" -j2 >"$WORK/bowtie2_make.log" 2>&1 || { tail -20 "$WORK/bowtie2_make.log"; echo "bowtie2 build failed"; }
  [ -x "$WORK/bowtie2/bowtie2-build" ] && BT2B="$WORK/bowtie2/bowtie2-build"
fi
[ -n "$BWA" ] || echo "bwa unavailable (set BWA_BIN) - skipped"
[ -n "$BT2B" ] || echo "bowtie2-build unavailable (set BT2_BUILD_BIN) - skipped"
mkdir -p "$WORK/idx"
BWAV=$( { "$BWA" 2>&1 || true; } | awk '/^Version:/{print $2}')
BTV=$( { "$BT2B" --version 2>&1 || true; } | awk '/version/{print $NF; exit}')
for x in ${INPUTS_ALN:-dna_1m dna_10m rep_10m dna_20m}; do
  n=$(awk -F'\t' -v k=$x '$1==k{print $2}' "$BENCH/MANIFEST.tsv")
  [ -z "$BWA" ] || meas --tool bwa-index --version "bwa $BWAV index -a bwtsw" --input $x --phase index --timeout "${TMO:-3600}" -- \
      python3 "$HERE/aligners/idx_wrap.py" "$WORK/idx/bwa_$x" "$n" -- "$BWA" index -a bwtsw -p "$WORK/idx/bwa_$x" "data/$x.fa"
  [ -z "$BT2B" ] || meas --tool bowtie2-build --version "bowtie2-build $BTV --threads 1" --input $x --phase index --timeout "${TMO:-3600}" -- \
      python3 "$HERE/aligners/idx_wrap.py" "$WORK/idx/bt2_$x" "$n" -- "$BT2B" --threads 1 --quiet "data/$x.fa" "$WORK/idx/bt2_$x"
done
echo "done: rows appended to $RAW"
