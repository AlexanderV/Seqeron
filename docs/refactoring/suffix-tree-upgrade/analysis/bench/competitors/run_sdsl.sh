#!/usr/bin/env bash
# sdsl-lite v2.1.1 (GPL-3.0): compressed suffix arrays / FM-indexes (csa_*) and compressed suffix tree (cst_sct3).
# Builds sdsl-lite from source (cmake, bundled submodules) then compiles sdsl_bench.cpp against it.
# Needs: git cmake g++ make
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
need git; need cmake; need g++; need make
SD_URL=https://github.com/simongog/sdsl-lite; SD_SHA=0546faf0552142f06ff4b201b671a5769dd007ad; SD_VER="sdsl-lite v2.1.1 (0546faf)"
clone_at $SD_URL $SD_SHA "$WORK/sdsl-lite"; note_licence "$WORK/sdsl-lite"
git -C "$WORK/sdsl-lite" submodule update --init --recursive --quiet || true
PREFIX="$WORK/sdsl-prefix"
if [ ! -f "$PREFIX/lib/libsdsl.a" ]; then
  # gcc>=13 no longer includes <cstdint> transitively
  ( cd "$WORK/sdsl-lite" && CXXFLAGS="-include cstdint" CFLAGS="-include stdint.h" ./install.sh "$PREFIX" >"$WORK/sdsl_install.log" 2>&1 ) || { tail -30 "$WORK/sdsl_install.log"; exit 1; }
fi
g++ -O3 -DNDEBUG -std=c++14 -include cstdint -I"$PREFIX/include" -L"$PREFIX/lib" "$HERE/sdsl/sdsl_bench.cpp" \
    -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -o "$WORK/sdsl_bench"
GV=$(g++ -dumpfullversion)
INPUTS_SDSL="${INPUTS_SDSL:-dna_1m dna_10m rep_10m prot_1m}"
TMO="${TMO:-2400}"
SKIP="${SKIP:-}"   # space separated cfg:input pairs to skip (already measured)
# The RRR-compressed FM-index and the CST are slow (locate = ~256 LF steps per hit at sample rate 512):
# for inputs >= 10 Mbp they run REPEATS_SLOW times (default 1) instead of REPEATS; the notes column records it via the repeats column.
REPEATS_SLOW="${REPEATS_SLOW:-1}"
for x in $INPUTS_SDSL; do
  for cfg in csa_wt_huff csa_sada csa_wt_rrr cst_sct3; do
    case " $SKIP " in *" $cfg:$x "*) continue;; esac
    extra=""; [ "$cfg" = cst_sct3 ] && [ "$x" = dna_1m ] && extra="lrs"
    r=$REPEATS
    case "$x" in dna_1m|prot_1m) ;; *) case $cfg in csa_wt_rrr|cst_sct3) r=$REPEATS_SLOW;; esac;; esac
    REPEATS=$r meas --tool "sdsl-$cfg" --version "$SD_VER g++-$GV -O3" --input $x --phase all --timeout "$TMO" -- "$WORK/sdsl_bench" $cfg data/$x.txt queries/$x $extra
  done
done
echo "done: rows appended to $RAW"
