#!/bin/bash
# E7 protocol (outside the shipped code): the end-to-end pipeline PH0-PH5 at several sizes, on the inputs prepared by run_e7_setup.sh.
# usage: run_e7_protocol.sh <work dir> <rows file> <group> <names...>      names: dna_50m dna_100m dna_200m dna_500m genome_50m genome_100m genome_200m ...
# groups: A        series A, fixed parameters: --scale-n 3000000000 --rule-n 100000000, i.e. every size runs with the values the scale rule gives at n = 10^8 (K_s 34,953, N_base 559,241, R 139,810, K_sort 34,953,
#                  run spill 2,185 B; PH3 K 33,333; PH4 R 266,667, 2 B2 tasks; PH5 block, window, slots); E7_RUNS runs (default 3), label <name>.e7.A
#         B        series B, ratio emulation: --scale-n 3000000000 (the E6 rule, f = n / 3x10^9); E7_RUNS runs (default 3), label <name>.e7.B
#         AB       both series, interleaved: for run 1..E7_RUNS, for every name, series B then series A (drift of the host shows in both)
#         threads  series B with 1 and 2 threads, E7_THREAD_RUNS runs each (default 1), labels <name>.e7.B.t1 / .t2
#         athreads series A with 1 and 2 threads (same)
# Every run: dropped page cache (cold), fresh process, DOTNET_GCHeapHardLimit = 96 MiB, 4 threads unless the group says otherwise, the independent sampler sample_e7.py (RssAnon, Dirty, scratch size, /proc/<pid>/io),
# the program's own rows (tagged E6 by the program, rewritten to E7 here), then (outside the timed region) the SHA-256 of the five outputs against the references' hashes of run_e7_setup.sh, the statistics lines,
# the scratch directory, and the outputs are deleted. A run whose rows already end in run_done is skipped (restart safe); an interrupted run's rows are purged and the run is repeated cold.
set -u
W=$1; ROWS=$2; GROUP=$3; shift 3
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000
DLL=${EMPROTO_DLL:-/tmp/stbench/bin_e7_run/SuffixTree.EmProto.dll}
HERE=$(cd "$(dirname "$0")" && pwd)
TIMEOUT=${E7_TIMEOUT:-2700}
cd "$W" || exit 1
printf 'R\tE7\tprotocol\thead\t%s\t\tgit rev-parse HEAD of the checkout when the group started; group %s; binary %s\n' "$(git -C /home/user/Seqeron rev-parse HEAD)" "$GROUP" "$DLL" >> "$ROWS"

purge() { # label run: drop the rows of an interrupted run
  awk -F'\t' -v l="$1" -v r="$2" '!($3==l && $4==r)' "$ROWS" > "$ROWS.tmp" && mv "$ROWS.tmp" "$ROWS"
}

compare() { # label run name outprefix: SHA-256 of every output against the reference hashes
  L=$1; R=$2; T=$3; O=$4
  declare -A GOT
  for K in pk sa lcp ovf cld; do ( sha256sum "$O.$K" 2>/dev/null | cut -d' ' -f1 > "$O.$K.sha" ) & done
  wait
  for K in pk sa lcp ovf cld; do
    GOT[$K]=$(cat "$O.$K.sha"); WANT=$(awk -v k=$K '$1==k {print $2}' "$T.refs.sha")
    if [ "$WANT" = "-" ]; then V=n/a; elif [ -n "${GOT[$K]}" ] && [ "${GOT[$K]}" = "$WANT" ]; then V=yes; else V=NO; fi   # n/a: no reference of that kind (the r500 cross-check has no child table)
    printf 'R\tE7\t%s\t%s\tcmp_%s\t%s\t\tSHA-256 %s of %s against the reference hash in %s.refs.sha\n' "$L" "$R" "$K" "$V" "${GOT[$K]:0:16}" "$(basename "$O.$K")" "$T" >> "$ROWS"
  done
  if ! grep -q '^stat ' "$T.refs.sha"; then V=n/a; elif [ "$(grep -E '^(nodecount|internal|lrsdepth|lrsrank)=' "$O.stats" 2>/dev/null | sort | sed 's/^/stat /')" = "$(grep '^stat ' "$T.refs.sha")" ]; then V=yes; else V=NO; fi
  printf 'R\tE7\t%s\t%s\tcmp_stats\t%s\t\tNODE_COUNT / internal / LRS_DEPTH / LRS_RANK against %s.refs.sha\n' "$L" "$R" "$V" "$T" >> "$ROWS"
  rm -f "$O".*.sha
}

run() { # input label run scratch-dir extra-args... : one cold run, then the comparison with the references
  T=$1; L=$2; R=$3; SCR=$4; shift 4
  if grep -q -P "^R\tE7\t$L\t$R\trun_done\t" "$ROWS"; then echo "$L run $R already done"; return 0; fi
  purge "$L" "$R"
  O=out.$L.$R; PIDF=pid.$L.$R
  rm -rf "$SCR" $O.* $PIDF 2>/dev/null
  FREE=$(df -B1G --output=avail /tmp | tail -1 | tr -d ' ')
  printf 'R\tE7\t%s\t%s\tdisk_free_at_start\t%s\tGB\tdf of /tmp before the run\n' "$L" "$R" "$FREE" >> "$ROWS"
  sync; echo 3 > /proc/sys/vm/drop_caches
  python3 -I "$HERE/sample_e7.py" $PIDF "$ROWS" "$L" "$R" "$SCR" &
  SP=$!
  sleep 0.2
  dotnet "$DLL" e6 --input "$T.txt" --out "$O" --scratch "$SCR" --case "$L" --run "$R" --rows "$ROWS.prog" "$@" > /dev/null &
  PID=$!
  echo $PID > $PIDF
  ( exec > /dev/null 2>&1; for ((w = 0; w < TIMEOUT; w++)); do kill -0 $PID 2>/dev/null || exit 0; sleep 1; done; kill -9 $PID 2>/dev/null ) &
  WD=$!
  wait $PID
  RC=$?
  wait $WD 2>/dev/null
  wait $SP
  [ -f "$ROWS.prog" ] && { sed 's/^R\tE6\t/R\tE7\t/' "$ROWS.prog" >> "$ROWS"; rm -f "$ROWS.prog"; }
  printf 'R\tE7\t%s\t%s\texit_code\t%s\t\tprogram exit code (0 = healthy, RssAnon cap 256 MB held, scratch gone; 137 = killed by the %s s watchdog)\n' "$L" "$R" "$RC" "$TIMEOUT" >> "$ROWS"
  compare "$L" "$R" "$T" "$O"
  printf 'R\tE7\t%s\t%s\tscratch_left\t%s\t\tentries of %s after the run (the directory must be gone)\n' "$L" "$R" "$( [ -e "$SCR" ] && ls -A "$SCR" | wc -l || echo 0)" "$SCR" >> "$ROWS"
  printf 'R\tE7\t%s\t%s\trun_done\tyes\t\tall rows of this run are written\n' "$L" "$R" >> "$ROWS"
  echo "$L run $R rc=$RC $(tail -n 600 "$ROWS" | awk -F'\t' -v l="$L" -v r="$R" '$3==l && $4==r && ($5=="total_s" || $5=="rssanon_peak") {printf "%s=%.4g ", $5, $6}')"
  rm -rf "$SCR" $O.* $PIDF
}

SCRATCH=scratch.d
case $GROUP in
AB) for i in $(seq 1 ${E7_RUNS:-3}); do for T in "$@"; do
      run $T "$T.e7.B" $i $SCRATCH --threads 4 --scale-n 3000000000
      run $T "$T.e7.A" $i $SCRATCH --threads 4 --scale-n 3000000000 --rule-n 100000000
    done; done ;;
A)  for T in "$@"; do for i in $(seq 1 ${E7_RUNS:-3}); do run $T "$T.e7.A" $i $SCRATCH --threads 4 --scale-n 3000000000 --rule-n 100000000; done; done ;;
B)  for T in "$@"; do for i in $(seq 1 ${E7_RUNS:-3}); do run $T "$T.e7.B" $i $SCRATCH --threads 4 --scale-n 3000000000; done; done ;;
threads)
    for T in "$@"; do for TH in 1 2; do for i in $(seq 1 ${E7_THREAD_RUNS:-1}); do run $T "$T.e7.B.t$TH" $i $SCRATCH --threads $TH --scale-n 3000000000; done; done; done ;;
athreads)
    for T in "$@"; do for TH in 1 2; do for i in $(seq 1 ${E7_THREAD_RUNS:-1}); do run $T "$T.e7.A.t$TH" $i $SCRATCH --threads $TH --scale-n 3000000000 --rule-n 100000000; done; done; done ;;
*) echo "unknown group $GROUP"; exit 2 ;;
esac
echo protocol-done
