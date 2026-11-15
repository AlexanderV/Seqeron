#!/bin/bash
# E7 part C: the cache sweep of one input (dna_100m or genome_100m) with run_e7_cache.sh, in the "lean" configuration (anonymous memory lowered so that a true 0.5x is possible at all):
#   DOTNET_GCHeapHardLimit = 48 MiB, --spill-kb 16 --tuple-kb 8 --budget-mb 64 --ph5-budget-mb 24 --b2-workers 2, 4 threads, --scale-n 3000000000 (series B rule), the .NET runtime pre-read outside the sub-group.
# usage: run_e7_cache_sweep.sh <work dir> <rows file> <name> <cell>...        cell = <cellname>:<limit MB or max>:<raise MB or ->:<watchdog s>:<runs>
#   e.g.   U:max:-:900:3   C250:250:-:900:3   D53:53:1024:900:3
set -u
W=$1; ROWS=$2; T=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd)
LEAN=(--threads 4 --scale-n 3000000000 --spill-kb 16 --tuple-kb 8 --budget-mb 64 --ph5-budget-mb 24 --b2-workers 2)
for C in "$@"; do
  IFS=: read NAME LIM RAISE WD RUNS <<< "$C"
  [ "$LIM" = max ] && LB=max || LB=$((LIM * 1048576))
  RZ=""; [ "$RAISE" != "-" ] && RZ=$((RAISE * 1048576))
  for i in $(seq 1 "$RUNS"); do
    E7_RAISE=$RZ E7_WATCHDOG=$WD E7_GC=0x3000000 "$HERE/run_e7_cache.sh" "$W" "$ROWS" "$T" "$T.e7.C.$NAME" "$i" "$LB" "${LEAN[@]}" 2>&1 | grep -v "^\s*$" | grep -E "run [0-9]+ rc="
  done
done
echo sweep-done
