#!/bin/bash
# E4 extra (report only, not a gate): the bounded mode's gather PLCP[SA[i]] reads the 4 B/char PLCP file at random; this sweeps the page cache that the build may use, with a cgroup-v1 memory limit
# (run_e4_cgroup.sh) at 100 Mbp: limit = 160 MiB (anonymous memory and the streamed files) + rho x 4n bytes (the PLCP file). The base mode is the control (sequential, no random reads).
# usage: run_e4_cgroup_sweep.sh <work dir with dna_100m.pk/.sa/.oracle.*> <rows file> [reps]    environment: RHOS (default "3 1.5 1 0.75 0.5 0.25"), CONTROL_RHOS (default "0.25 0"), TIMEOUT (s per run, default 240):
# a run above the limit is killed and recorded as status TIMEOUT (random 4 KB reads of the PLCP file without the cache: hours).
set -u
W=$1; ROWS=$2; REPS=${3:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
export DOTNET_GCHeapHardLimit=0x6000000
cd "$W" || exit 1
mkdir -p scratch
N=$(grep '^n=' dna_100m.pk.meta | cut -d= -f2)
F=$((4 * N)); A=$((160 * 1048576))
R=266667
run() { # rho mode-args label rep
  RHO=$1; MODE=$2; LABEL=$3; REP=$4
  LIM=$(python3 -I -c "print(int($A + $RHO * $F))")
  timeout -k 5 "${TIMEOUT:-240}" bash "$HERE/run_e4_cgroup.sh" "$LIM" "$LABEL-$REP" --text dna_100m.pk --sa dna_100m.sa --out dna_100m.out --scratch scratch --scale-n 3000000000 --verify dna_100m.oracle $MODE --case "$LABEL" --run "$REP" --rows "$ROWS" 2>>"$ROWS.cgroup.log" > /dev/null
  RC=$?
  if [ $RC -eq 124 ] || [ $RC -eq 137 ]; then
    pkill -f "emproto_e4_$LABEL-$REP" ; pkill -f 'EmProto.dll e4 --cgroup-usage'; sleep 3; rmdir "/sys/fs/cgroup/memory/emproto_e4_$LABEL-$REP" 2>/dev/null; rm -rf scratch/*
    printf 'R\tE4\t%s\t%s\tstatus\tTIMEOUT\t\tkilled after %s s (limit %s bytes = 160 MiB + %s x 4n)\n' "$LABEL" "$REP" "${TIMEOUT:-240}" "$LIM" "$RHO" >> "$ROWS"
  fi
  echo "$LABEL rep $REP rc=$RC"
}
for rep in $(seq 1 "$REPS"); do
  for rho in ${RHOS:-3 1.5 1 0.75 0.5 0.25}; do run $rho "--mode bounded --rounds 1" "cg.b1.rho$rho" $rep; done
  for rho in ${CONTROL_RHOS:-0.25 0}; do run $rho "--mode distributed" "cg.dist.rho$rho" $rep; done
done
echo sweep-done
