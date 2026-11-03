#!/usr/bin/env bash
# libsais (Apache-2.0) and libdivsufsort (MIT): build from source with cmake, link sa_bench.c, measure.
# Needs: git cmake gcc (OpenMP optional but used for the libsais_omp rows)
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
need git; need cmake; need gcc
LS_URL=https://github.com/IlyaGrebnov/libsais;      LS_SHA=1aae6344b8a58f00626f001b1fe33bdbecd5e469;   LS_VER="libsais v2.9.1 (1aae634)"
DS_URL=https://github.com/y-256/libdivsufsort;     DS_SHA=82b01f764e342f1858d4473c17b56844aac9c12d;   DS_VER="libdivsufsort 2.0.1 (82b01f7)"
clone_at $LS_URL $LS_SHA "$WORK/libsais";           note_licence "$WORK/libsais"
clone_at $DS_URL $DS_SHA "$WORK/libdivsufsort";     note_licence "$WORK/libdivsufsort"

cmake -S "$WORK/libsais" -B "$WORK/libsais/build" -DCMAKE_BUILD_TYPE=Release -DLIBSAIS_USE_OPENMP=ON -DLIBSAIS_BUILD_SHARED_LIB=OFF >/dev/null
cmake --build "$WORK/libsais/build" -j2 >/dev/null
cmake -S "$WORK/libdivsufsort" -B "$WORK/libdivsufsort/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DBUILD_EXAMPLES=OFF \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 >/dev/null
cmake --build "$WORK/libdivsufsort/build" -j2 >/dev/null

LSLIB=$(find "$WORK/libsais/build" -name 'liblibsais.a' -o -name 'libsais.a' | head -1)
DSLIB=$(find "$WORK/libdivsufsort/build" -name 'libdivsufsort.a' | head -1)
gcc -O3 -fopenmp -DUSE_LIBSAIS -DLIBSAIS_OPENMP -DLIBSAIS_OMP -I"$WORK/libsais/include" "$HERE/sa_c/sa_bench.c" "$LSLIB" -lpthread -o "$WORK/sa_bench_libsais"
gcc -O3 -DUSE_DIVSUFSORT -I"$WORK/libdivsufsort/build/include" "$HERE/sa_c/sa_bench.c" "$DSLIB" -lpthread -o "$WORK/sa_bench_divsufsort"

GCCV=$(gcc -dumpfullversion)
INPUTS_SA="${INPUTS_SA:-$INPUTS_DNA dna_50m prot_1m}"
for x in $INPUTS_SA; do
  meas --tool libsais        --version "$LS_VER gcc-$GCCV -O3 (OpenMP build; serial rows use libsais())" --input $x --phase all -- "$WORK/sa_bench_libsais" data/$x.txt queries/$x
  meas --tool libdivsufsort  --version "$DS_VER gcc-$GCCV -O3"                                           --input $x --phase all -- "$WORK/sa_bench_divsufsort" data/$x.txt queries/$x
done
echo "done: rows appended to $RAW ; aggregate with: python3 $BENCH/lib/assemble.py $BENCH"
