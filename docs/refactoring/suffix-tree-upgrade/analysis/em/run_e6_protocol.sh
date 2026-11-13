#!/bin/bash
# E6 protocol (outside the shipped code): the end-to-end pipeline PH0-PH5 on the inputs prepared by run_e6_setup.sh.
# usage: run_e6_protocol.sh <work dir with <name>.txt and the references> <rows file> <group> [names...]
# groups: main    (G1/G2/G3: 4 threads, cold cache, scale-down by n / 3e9, 3 runs per input; every run is compared byte for byte with the references by cmp)
#         threads (1 and 2 threads, 3 runs each, cmp against the same references: identity over thread counts and the scaling)
#         kvar    (two other values of K: K/2 and 2K for PH3 and PH1, 4 threads, 1 run each, cmp)
#         dirty   (the scratch on tmpfs, so that the system-wide Dirty of /proc/meminfo is the dirty pages of the outputs only; 3 runs)
#         small   (dna_10m rep_10m prot_1m and the families: cmp against the references, SA hash against the committed libsais hash where there is one)
# Every run: dropped page cache (cold), fresh process, DOTNET_GCHeapHardLimit = 96 MiB, real block sizes, the independent sampler sample_e6.py (RssAnon, Dirty, scratch size, /proc/<pid>/io).
set -u
W=$1; ROWS=$2; GROUP=$3; shift 3
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e6_run/SuffixTree.EmProto.dll}
HERE=$(cd "$(dirname "$0")" && pwd)
VECTORS=$HERE/../sa_vectors.tsv
cd "$W" || exit 1
printf 'R\tE6\tprotocol\thead\t%s\t\tgit rev-parse HEAD of the binary under test; group %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$GROUP" >> "$ROWS"

cmp_row() { # label run what produced reference
  if cmp -s "$4" "$5"; then V=yes; else V=NO; fi
  printf 'R\tE6\t%s\t%s\tcmp_%s\t%s\t\tcmp of %s against %s\n' "$1" "$2" "$3" "$V" "$(basename "$4")" "$(basename "$5")" >> "$ROWS"
}

compare() { # label run name outprefix: cmp every output against the references of the setup
  L=$1; R=$2; T=$3; O=$4
  cmp_row "$L" "$R" pk "$O.pk" "$T.ref.pk"
  cmp_row "$L" "$R" sa "$O.sa" "$T.refsa"
  cmp_row "$L" "$R" lcp "$O.lcp" "$T.oracle.lcp"
  cmp_row "$L" "$R" ovf "$O.ovf" "$T.oracle.ovf"
  cmp_row "$L" "$R" cld "$O.cld" "$T.ref.cldP"
  if diff <(grep -E '^(nodecount|internal|lrsdepth|lrsrank)=' "$O.stats") <(grep -E '^(nodecount|internal|lrsdepth|lrsrank)=' "$T.ref.stats") > /dev/null; then V=yes; else V=NO; fi
  printf 'R\tE6\t%s\t%s\tcmp_stats\t%s\t\tNODE_COUNT / internal / LRS_DEPTH / LRS_RANK against %s.ref.stats\n' "$L" "$R" "$V" "$T" >> "$ROWS"
}

run() { # input label run scratch-dir extra-args... : one cold run, then the comparison with the references
  T=$1; L=$2; R=$3; SCR=$4; shift 4
  O=out.$T.$L.$R
  rm -rf "$SCR" $O.* $T.pidfile 2>/dev/null
  mkdir -p "$(dirname "$SCR")"
  sync; echo 3 > /proc/sys/vm/drop_caches
  python3 -I "$HERE/sample_e6.py" $T.pidfile "$ROWS" "$L" "$R" "$SCR" &
  SP=$!
  sleep 0.2
  ( echo $BASHPID > $T.pidfile; exec dotnet "$DLL" e6 --input "$T.txt" --out "$O" --scratch "$SCR" --case "$L" --run "$R" --rows "$ROWS" "$@" > /dev/null )
  RC=$?
  wait $SP
  printf 'R\tE6\t%s\t%s\texit_code\t%s\t\tprogram exit code (0 = healthy, RssAnon cap held, scratch gone, every --verify-* equal)\n' "$L" "$R" "$RC" >> "$ROWS"
  compare "$L" "$R" "$T" "$O"
  if [ -n "${VEC:-}" ]; then
    H=$(DOTNET_GCHeapHardLimit= dotnet "$DLL" hash-sa --file "$O.sa" --case "$L" | cut -f5)
    E=$(awk -F'\t' -v c="$VEC" '$1==c {print $5}' "$VECTORS")
    if [ "$H" = "$E" ]; then V=yes; else V=NO; fi
    printf 'R\tE6\t%s\t%s\tsa_equals_libsais\t%s\t\tSHA-256 of SA[1..n] %s against the committed libsais hash of corpus %s (analysis/sa_vectors.tsv)\n' "$L" "$R" "$V" "$H" "$VEC" >> "$ROWS"
  fi
  printf 'R\tE6\t%s\t%s\tscratch_left\t%s\t\tentries of %s after the run (the directory must be gone)\n' "$L" "$R" "$( [ -e "$SCR" ] && ls -A "$SCR" | wc -l || echo 0)" "$SCR" >> "$ROWS"
  echo "$L run $R rc=$RC"
  rm -f $O.* $T.pidfile
}

BASE_ARGS=(--scale-n 3000000000)
set_names() { NAMES=("$@"); }
case $GROUP in
main)
  [ $# -gt 0 ] && set_names "$@" || set_names dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    for i in 1 2 3; do
      if [ $i = 1 ]; then V=(--verify-pk $T.ref.pk --verify-sa $T.refsa --verify-lcp $T.oracle --verify-child $T.ref.cldP --verify-stats $T.ref.stats); else V=(); fi
      run $T "$T.e6.t4" $i scratch.d --threads 4 "${BASE_ARGS[@]}" "${V[@]}"
    done
  done ;;
threads)
  [ $# -gt 0 ] && set_names "$@" || set_names dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    for TH in 1 2; do
      for i in 1 2 3; do run $T "$T.e6.t$TH" $i scratch.d --threads $TH "${BASE_ARGS[@]}"; done
    done
  done ;;
kvar)
  [ $# -gt 0 ] && set_names "$@" || set_names dna_100m genome_100m
  for T in "${NAMES[@]}"; do
    # the scale-down values are K = 33,333 (PH3) and K_s = 34,953 (PH1); half and double
    run $T "$T.e6.t4.Khalf" 1 scratch.d --threads 4 "${BASE_ARGS[@]}" --k 16666 --k-s 17476
    run $T "$T.e6.t4.Kdouble" 1 scratch.d --threads 4 "${BASE_ARGS[@]}" --k 66666 --k-s 69905
  done ;;
dirty)
  [ $# -gt 0 ] && set_names "$@" || set_names dna_100m genome_100m
  mkdir -p /dev/shm/e6
  for T in "${NAMES[@]}"; do
    for i in 1 2 3; do run $T "$T.e6.t4.shm" $i /dev/shm/e6/scratch.d --threads 4 "${BASE_ARGS[@]}"; done
  done ;;
small)
  [ $# -gt 0 ] && set_names "$@" || set_names dna_10m rep_10m prot_1m fam_a_1e5 fam_fib_1e5 fam_tm_1e5 fam_tandem_1e5 fam_a_1e6 fam_fib_1e6 fam_tm_1e6 fam_nrun_1e6 fam_genome_1e6
  for T in "${NAMES[@]}"; do
    case $T in
      dna_10m|rep_10m|prot_1m) VEC=$T ;;
      fam_a_1e5) VEC=all-same-n100000 ;; fam_a_1e6) VEC=all-same-n1000000 ;;
      fam_fib_1e5) VEC=fibonacci-n100000 ;; fam_fib_1e6) VEC=fibonacci-n1000000 ;;
      fam_tm_1e5) VEC=thue-morse-n100000 ;; fam_tm_1e6) VEC=thue-morse-n1000000 ;;
      *) VEC= ;;
    esac
    run $T "$T.e6.t4" 1 scratch.d --threads 4 "${BASE_ARGS[@]}"
  done
  VEC= ;;
*) echo "unknown group $GROUP"; exit 2 ;;
esac
echo protocol-done
