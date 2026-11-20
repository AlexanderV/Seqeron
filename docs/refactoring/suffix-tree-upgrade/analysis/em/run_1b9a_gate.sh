#!/bin/bash
# Gate of work package 1b.9a (outside the shipped code): the file-backed tree (Checked, and Trusted for the record) against the in-memory tree, count of present 20-mers, locate of 8-mers and count of absent 24-mers, one tree per
# process, warm cache, the variants in a rotated order from round to round (run_1b9a_matrix.sh), medians over the processes (summarize_1b9a.py). The machine must be idle; every step builds nothing but the containers.
#   usage: run_1b9a_gate.sh <work dir> <step>...   steps:
#     harness                 build the harness (Release) into <work>/bin/new
#     containers              build the containers (and a small one of 1 Mbp for the JIT warm-up of the cold open) of the texts below with the builder of the head (the harness variant `container` of SesaBuildCostTests; the EmProto.Tests project must be built)
#     m10 | m5 | mprot | m100 the matrix at 10 Mbp of DNA / 10 Mbp over five symbols / 10 Mbp protein-like / 100 Mbp of DNA (ROUNDS_10 default 15, ROUNDS_100 default 7, memory and Checked only at 100 Mbp)
#     open                    the open cost: warm (25 opens in a JIT-warm process) and cold (page cache dropped before each of 10 processes; a JIT-warm process, the small container warms the code)
#     inmem                   the in-memory tree of the head against the one of another build (OLD_BIN, e4a0b2da built with -p:NoOpen=true), 15 rounds, memory only
#   texts: DNA_10M, RND5_10M, PROT_10M, DNA_100M (paths in the environment, defaults below); the 100 Mbp container is built with the scaled geometry (G3_SCALE 3000000000), the others with the real one.
set -u
W=$(readlink -f -m "$1"); shift
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
DNA_10M=${DNA_10M:-/tmp/stbench/data/dna_10m.txt}
RND5_10M=${RND5_10M:-/tmp/stbench/g9a/rnd5_10m.txt}
PROT_10M=${PROT_10M:-/tmp/stbench/g9a/prot_10m.txt}
DNA_100M=${DNA_100M:-/tmp/stbench/e6/dna_100m.txt}
ROUNDS_10=${ROUNDS_10:-15}
ROUNDS_100=${ROUNDS_100:-7}
mkdir -p "$W/rows" "$W/bin"
build_container() { # text target scale
  rm -f "$3.tmp"; rm -f "$2"
  env DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1 G3_VARIANT=container G3_INPUT=$1 G3_TARGET=$2 G3_ROWS=$W/rows/build.tsv G3_LABEL="$(basename "$2")" G3_SCALE=$3 G3_PSI=1 \
    dotnet test "$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaBuildCostTests" > "$W/build.log" 2>&1 || { echo "container build failed: $2"; tail -20 "$W/build.log"; exit 1; }
}
matrix() { # label text container alphabet rounds [variants]
  "$HERE/run_1b9a_matrix.sh" "$W/rows/$1.tsv" "$5" "${6:-memory checked trusted}" "$W/bin/new" "$3" "$2" "$4"
  python3 "$HERE/summarize_1b9a.py" "$W/rows/$1.tsv" | tee "$W/rows/$1.summary.txt"
}
for STEP in "$@"; do
  case $STEP in
    harness) dotnet build "$ROOT/docs/refactoring/suffix-tree-upgrade/analysis/em/readbench/ReadBench.csproj" -c Release -o "$W/bin/new" 2>&1 | grep -E "error|rror\(s\)" ;;
    containers)
      build_container "$DNA_10M" "$W/dna_10m.sesa" 0
      build_container "$RND5_10M" "$W/rnd5_10m.sesa" 0
      build_container "$PROT_10M" "$W/prot_10m.sesa" 0
      build_container "${DNA_1M:-/tmp/stbench/data/dna_1m.txt}" "$W/small.sesa" 0
      ls -l "$W"/*.sesa ;;
    m10) matrix dna_10m "$DNA_10M" "$W/dna_10m.sesa" ACGT "$ROUNDS_10" ;;
    m5) matrix rnd5_10m "$RND5_10M" "$W/rnd5_10m.sesa" ACGTN "$ROUNDS_10" ;;
    mprot) matrix prot_10m "$PROT_10M" "$W/prot_10m.sesa" ACDEFGHIKLMNPQRSTVWY "$ROUNDS_10" ;;
    m100)
      build_container "$DNA_100M" "$W/dna_100m.sesa" 3000000000
      # the new codec at scale: the independent checker of e8-verify (eleven checks over the sections read in place) and the Structure level of the reader, before the throughput
      dotnet "$ROOT/apps/SuffixTree.EmProto/bin/Release/net10.0/SuffixTree.EmProto.dll" e8-verify --container "$W/dna_100m.sesa" --input "$DNA_100M" --threads 4 --scratch "$W" --case dna_100m.verify --rows "$W/rows/verify_100m.tsv" > "$W/rows/verify_100m.log" 2>&1
      printf '#e8-verify\tdna_100m\texit\t%s\n' "$?" >> "$W/rows/verify_100m.tsv"
      env -u DOTNET_GCHeapHardLimit SESA_ROWS="$W/rows/structure_100m.tsv" SESA_LABEL=dna_100m.structure SESA_VARIANT=open SESA_LEVEL=Structure SESA_ACCESS=Checked SESA_FILE="$W/dna_100m.sesa" SESA_INPUT="$DNA_100M" \
        dotnet test "$ROOT/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaReaderCostTests" > "$W/rows/structure_100m.log" 2>&1
      matrix dna_100m "$DNA_100M" "$W/dna_100m.sesa" ACGT "$ROUNDS_100" "memory checked"
      rm -f "$W/dna_100m.sesa" ;;
    open)
      for i in 1 2 3; do dotnet "$W/bin/new/SuffixTree.Esa.Bench.dll" open "$W/dna_10m.sesa" x; done | tee "$W/rows/open_warm.tsv"
      : > "$W/rows/open_cold.tsv"
      for i in $(seq 1 10); do sync; echo 3 > /proc/sys/vm/drop_caches; dotnet "$W/bin/new/SuffixTree.Esa.Bench.dll" opencold "$W/dna_10m.sesa" "$W/small.sesa" >> "$W/rows/open_cold.tsv"; done
      cat "$W/rows/open_cold.tsv" ;;
    inmem)
      OLD=${OLD_BIN:?OLD_BIN: the harness build of the other commit (e4a0b2da, built with -p:NoOpen=true)}
      "$HERE/run_1b9a_matrix.sh" "$W/rows/inmem_old_new.tsv" 15 "memory@$OLD memory@$W/bin/new" "$W/bin/new" "$W/dna_10m.sesa" "$DNA_10M" ACGT
      python3 "$HERE/summarize_1b9a.py" "$W/rows/inmem_old_new.tsv" memory@old | tee "$W/rows/inmem_old_new.summary.txt" ;;
    *) echo "unknown step $STEP"; exit 2 ;;
  esac
done
