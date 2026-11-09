#!/bin/bash
# E3 protocol (outside the shipped code): baseline, thread/K determinism, and the page-cache sweep under a cgroup-v1 memory limit.
# usage: run_e3_protocol.sh <label> <packed text> <rank file> <oracle sa hash> <rows file> [reps]
# environment: RHOS="2.0 1.0 0.5 0.25 0" (cache ratios to sweep), ONLY_RHO=1 (skip the baseline and the determinism runs)
# limit = A + rho * (text + rank) + S  with A = 85 MB (anonymous allowance: above the 83 MB peak RssAnon of the unlimited run), S = 8 MB (streaming allowance)
set -u
LABEL=$1; TEXT=$2; RANKS=$3; ORACLE=$4; ROWS=$5; REPS=${6:-3}
HERE=$(dirname "$(readlink -f "$0")")
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_GCHeapHardLimit=0x6000000          # 96 MiB managed heap cap (ADR-004 section 5)
export EMPROTO_DLL=${EMPROTO_DLL:-/home/user/Seqeron/apps/SuffixTree.EmProto/bin/Release/net10.0/SuffixTree.EmProto.dll}
SA=${SA:-/tmp/stbench/e3/sa.bin}; SCR=${SCR:-/tmp/stbench/e3/scratch}
MB=1048576
T=$(( $(stat -c %s "$TEXT") + $(stat -c %s "$RANKS") ))
verify() { # case
  H=$(dotnet "$EMPROTO_DLL" hash-sa --file "$SA" --case "$1" --rows "$ROWS" | cut -f5)
  if [ "$H" = "$ORACLE" ]; then R=yes; else R=NO; fi
  printf 'R\tE3\t%s\t%s\tsa_equals_oracle\t%s\t\tSHA-256 of SA[1..n] against SaisKernel in RAM on the same text\n' "$1" "$2" "$R" | tee -a "$ROWS" > /dev/null
  echo "$1 run $2 sa_equals_oracle=$R"
}
run() { # limit case run extra-args...
  LIM=$1; CASE=$2; RUN=$3; shift 3
  rm -f "$SA"
  bash "$HERE/run_e3_cgroup.sh" "$LIM" "$CASE-$RUN" --text "$TEXT" --ranks "$RANKS" --sa "$SA" --scratch "$SCR" --case "$CASE" --run "$RUN" --rows "$ROWS" "$@" 2>>"$ROWS.cgroup.log" >/dev/null
  RC=$?
  if [ $RC -ne 0 ]; then
    # the build did not finish (rc 137 = killed by the OOM killer of the sub-group): its SA file is incomplete, there is nothing to compare
    printf 'R\tE3\t%s\t%s\tstatus\tFAILED\t\tprocess exit code %s (137 = SIGKILL by the cgroup OOM killer); no hash comparison\n' "$CASE" "$RUN" "$RC" | tee -a "$ROWS" > /dev/null
    echo "$CASE run $RUN FAILED rc=$RC"
    return
  fi
  verify "$CASE" "$RUN"
}
if [ -z "${ONLY_RHO:-}" ]; then
  for i in $(seq 1 "$REPS"); do run max "$LABEL.baseline_4t_K3000" "$i" --threads 4; done
  for th in 1 2 4; do run max "$LABEL.threads$th" 1 --threads "$th"; done
  run max "$LABEL.K1500_4t" 1 --threads 4 --k-div 1500
  run max "$LABEL.K6000_4t" 1 --threads 4 --k-div 6000
fi
for rho in ${RHOS:-2.0 1.0 0.5 0.25 0}; do
  LIM=$(python3 -I -c "print(int((85 + 8) * $MB + $rho * $T))")
  for i in $(seq 1 "$REPS"); do run "$LIM" "$LABEL.cache_rho$rho" "$i" --threads 4; done
done
echo protocol-done
