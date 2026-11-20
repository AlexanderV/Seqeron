#!/bin/bash
# Benchmark definition of work package 1b.9b (BENCHMARK-1b9b.md in analysis/em/readbench): the in-memory tree against the file-backed tree (Checked) on six texts of 10 Mbp and DNA at 100 Mbp, with two instruments, always both.
#   (a) cross-process: one tree per process, the variants in a rotated order from round to round, ROUNDS processes per variant (default 21), warm page cache, no heap limit (run_1b9a_matrix.sh, summarized by summarize_1b9b.py);
#   (b) in process: the trees side by side in one process, the passes in a rotated order (ReadBench variant memory+checked+trusted), PROCS processes of PASSES passes.
# The machine must be idle: no build, no test, no other timing run overlaps.
#   usage: run_1b9b_bench.sh <work dir> <step>...   steps:
#     harness              build the harness (Release) into <work>/bin/head
#     texts                generate the five texts that have a generator into <work>/data (the DNA text is the S0 input DNA_10M)
#     containers           build the containers of the six texts (and a small one) with the builder of the head (SesaBuildCostTests, variant container; the EmProto.Tests project must be built in Release)
#     cross  <text>...     the cross-process instrument on these texts (dna sigma5 protu protr text90 bytes256), rows to <work>/rows/cross_<text>.tsv
#     inproc <text>...     the in-process instrument on these texts, rows to <work>/rows/inproc_<text>.tsv
#     report               the table of both instruments over every rows file of <work>/rows
#     m100                 DNA at 100 Mbp (DNA_100M, default /tmp/stbench/e6/dna_100m.txt): the container with the scaled geometry of the builder (G3_SCALE 3000000000), the cross-process instrument (ROUNDS_100 rounds, default 7) and the in-process instrument (IP100 processes, default 4, 5 passes); the container is removed unless KEEP100 is set
#   isolated experiments (1b.9b; each needs the harness built with -p:Experiments=true into <work>/bin/headx for the ones marked X):
#     null <text>...       the null experiment of the cross-process instrument: the in-memory tree against itself (a copy of the harness directory, <work>/bin/head2): the ratio of the medians of two identical variants is the noise floor
#     jit <text>...        the cross-process in-memory and file-backed trees under three JIT settings (default, TieredPGO=0, TieredCompilation=0): is the spread between processes the tiering?
#     abl <text>...    X   the in-process ablation at the engine level: the facade, the index, the specialised descent over pinned arrays (byte text, packed text copy) of the in-memory tree
#     old <text>...        the in-memory tree of the head against the one of another commit (OLD_DLL: the SuffixTree.Esa.dll of e4a0b2da, built with its own dependencies) and against the head before the change (PRE_DLL), in process
#     threads <text>...  X the four-thread instrument: one tree, four threads, the facade against the index (the reader count), memory and file-backed, in process (RB_THREADS, RB_REPEAT; ROUNDS rounds)
#   environment: ROUNDS (21), PROCS (6), PASSES (11), DNA_10M, RB_L20 / RB_L8 / RB_L24 (pattern lengths 20 / 8 / 24), RB_N20 / RB_N8 / RB_N24 (numbers 200000 / 20000 / 100000), TAG (a suffix of the rows file names)
set -u
W=$(readlink -f -m "$1"); shift
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
DNA_10M=${DNA_10M:-/tmp/stbench/data/dna_10m.txt}
ROUNDS=${ROUNDS:-21}; PROCS=${PROCS:-6}; PASSES=${PASSES:-11}; TAG=${TAG:-}
BIN=${BIN:-$W/bin/head}
mkdir -p "$W/rows" "$W/bin" "$W/data"
# text name -> file, alphabet argument of the harness
text_file() { case $1 in dna) echo "$DNA_10M";; *) echo "$W/data/$1_10m.txt";; esac; }
text_alpha() { case $1 in dna) echo ACGT;; sigma5) echo ACGTN;; protu|protr) echo ACDEFGHIKLMNPQRSTVWY;; *) echo @text;; esac; }
build_container() { # text target scale
  rm -f "$2" "$2.tmp"
  env DOTNET_GCHeapHardLimit=0x6000000 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1 G3_VARIANT=container G3_INPUT=$1 G3_TARGET=$2 G3_ROWS=$W/rows/build.tsv G3_LABEL="$(basename "$2")" G3_SCALE=$3 G3_PSI=1 \
    dotnet test "$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj" -c Release --no-build --nologo --filter "FullyQualifiedName~SesaBuildCostTests" > "$W/build.log" 2>&1 || { echo "container build failed: $2"; tail -20 "$W/build.log"; exit 1; }
}
run_texts() { # kind text...
  kind=$1; shift
  for t in "$@"; do
    file=$(text_file "$t"); alpha=$(text_alpha "$t"); rows="$W/rows/${kind}_${t}${TAG}.tsv"
    if [ "$kind" = cross ]; then
      "$HERE/run_1b9a_matrix.sh" "$rows" "$ROUNDS" "memory checked" "$BIN" "$W/${t}_10m.sesa" "$file" "$alpha"
    else
      : > "$rows"
      for p in $(seq 1 "$PROCS"); do sync; dotnet "$BIN/SuffixTree.Esa.Bench.dll" memory+checked+trusted "$W/${t}_10m.sesa" "$file" "$PASSES" "ip$p" "$alpha" >> "$rows"; done
    fi
  done
}
while [ $# -gt 0 ]; do
  STEP=$1; shift
  case $STEP in
    harness) dotnet build "$ROOT/docs/refactoring/suffix-tree-upgrade/analysis/em/readbench/ReadBench.csproj" -c Release -o "$W/bin/head" 2>&1 | grep -E " error |rror\(s\)" ;;
    texts)
      for spec in "sigma5 sigma5" "prot-uniform protu" "prot-real protr" "text90 text90" "bytes256 bytes256"; do
        set -- $spec "$@"; k=$1; n=$2; shift 2
        dotnet "$BIN/SuffixTree.Esa.Bench.dll" gen "$k" 10000000 20261120 "$W/data/${n}_10m.txt"
      done ;;
    containers)
      for t in dna sigma5 protu protr text90 bytes256; do build_container "$(text_file $t)" "$W/${t}_10m.sesa" 0; done
      build_container /tmp/stbench/data/dna_1m.txt "$W/small.sesa" 0
      ls -l "$W"/*.sesa ;;
    cross|inproc)
      names=()
      while [ $# -gt 0 ]; do case $1 in harness|texts|containers|cross|inproc|report|null|jit|abl|old|threads|m100) break;; *) names+=("$1"); shift;; esac; done
      run_texts "$STEP" "${names[@]}" ;;
    null)
      names=()
      while [ $# -gt 0 ]; do case $1 in harness|texts|containers|cross|inproc|report|null|jit|abl|old|threads|m100) break;; *) names+=("$1"); shift;; esac; done
      [ -d "$W/bin/head2" ] || cp -r "$W/bin/head" "$W/bin/head2"
      for t in "${names[@]}"; do
        "$HERE/run_1b9a_matrix.sh" "$W/rows/null_${t}${TAG}.tsv" "$ROUNDS" "memory memory@$W/bin/head2" "$BIN" "$W/${t}_10m.sesa" "$(text_file $t)" "$(text_alpha $t)"
      done ;;
    jit)
      names=()
      while [ $# -gt 0 ]; do case $1 in harness|texts|containers|cross|inproc|report|null|jit|abl|old|threads|m100) break;; *) names+=("$1"); shift;; esac; done
      for t in "${names[@]}"; do
        for setting in default DOTNET_TieredPGO:0 DOTNET_TieredCompilation:0; do
          env_spec=""; [ "$setting" = default ] || env_spec="=$setting"
          "$HERE/run_1b9a_matrix.sh" "$W/rows/jit_${t}_${setting//[:=]/_}${TAG}.tsv" "${JIT_ROUNDS:-12}" "memory$env_spec checked$env_spec" "$BIN" "$W/${t}_10m.sesa" "$(text_file $t)" "$(text_alpha $t)"
        done
      done ;;
    abl)
      names=()
      while [ $# -gt 0 ]; do case $1 in harness|texts|containers|cross|inproc|report|null|jit|abl|old|threads|m100) break;; *) names+=("$1"); shift;; esac; done
      for t in "${names[@]}"; do
        rows="$W/rows/abl_${t}${TAG}.tsv"; : > "$rows"
        for p in $(seq 1 "$PROCS"); do sync; dotnet "$W/bin/headx/SuffixTree.Esa.Bench.dll" "${ABL_VARIANTS:-memory+ix-memory+ix-mem-b+ix-mem-p}" "$W/${t}_10m.sesa" "$(text_file $t)" "$PASSES" "ab$p" "$(text_alpha $t)" >> "$rows"; done
      done ;;
    old)
      names=()
      while [ $# -gt 0 ]; do case $1 in harness|texts|containers|cross|inproc|report|null|jit|abl|old|threads|m100) break;; *) names+=("$1"); shift;; esac; done
      for t in "${names[@]}"; do
        rows="$W/rows/old_${t}${TAG}.tsv"; : > "$rows"
        for p in $(seq 1 "$PROCS"); do sync; dotnet "$BIN/SuffixTree.Esa.Bench.dll" "alc=${OLD_DLL:?OLD_DLL}+memory+alc=${PRE_DLL:?PRE_DLL}" "$W/${t}_10m.sesa" "$(text_file $t)" "$PASSES" "o$p" "$(text_alpha $t)" >> "$rows"; done
      done ;;
    threads)
      names=()
      while [ $# -gt 0 ]; do case $1 in harness|texts|containers|cross|inproc|report|null|jit|abl|old|threads|m100) break;; *) names+=("$1"); shift;; esac; done
      for t in "${names[@]}"; do
        rows="$W/rows/threads_${t}${TAG}.tsv"; : > "$rows"
        for p in $(seq 1 "${TPROCS:-3}"); do sync; dotnet "$W/bin/headx/SuffixTree.Esa.Bench.dll" "mt:${MT_VARIANTS:-memory+ix-memory+checked+ix-file}" "$W/${t}_10m.sesa" "$(text_file $t)" "$ROUNDS" "mt$p" "$(text_alpha $t)" >> "$rows"; done
      done ;;
    m100)
      DNA_100M=${DNA_100M:-/tmp/stbench/e6/dna_100m.txt}
      build_container "$DNA_100M" "$W/dna_100m.sesa" 3000000000
      "$HERE/run_1b9a_matrix.sh" "$W/rows/cross_dna100${TAG}.tsv" "${ROUNDS_100:-7}" "memory checked" "$BIN" "$W/dna_100m.sesa" "$DNA_100M" ACGT
      : > "$W/rows/inproc_dna100${TAG}.tsv"
      for p in $(seq 1 "${IP100:-4}"); do sync; dotnet "$BIN/SuffixTree.Esa.Bench.dll" memory+checked+trusted "$W/dna_100m.sesa" "$DNA_100M" 5 "ip$p" ACGT >> "$W/rows/inproc_dna100${TAG}.tsv"; done
      [ -n "${KEEP100:-}" ] || rm -f "$W/dna_100m.sesa" ;;
    report) python3 "$HERE/summarize_1b9b.py" "$W/rows" ;;
    *) echo "unknown step $STEP"; exit 2 ;;
  esac
done
