#!/usr/bin/env bash
# The A/B session of work package S2.3 (ADR-005 section 3.4; the instrument is analysis/bench/abapp): the S1 head built to BASE (dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o $BASE)
# against the S2.3 head built to NEW, on the containers that the BASELINE builder wrote (the new builder writes the same bytes: results_s23_identity.txt), count (present and absent) and locate. One cell = one input, one kind of spec and one
# thread count; one cell = 5 fresh processes of abapp (5 ABBA blocks of 1.0 s per operation each); the figure of a process is the median of its per-block ratios B / A; the figure of a cell is the median of the 5 (summarize_s23_ab.py).
# The null experiment (A = B = BASE) uses the same cells. The libsais drift control (drift.sh of analysis/s2) brackets the session.
#   usage: run_s23_ab.sh <null|ab> <cell>...   with <cell> = <input>:<kind>:<threads>   (kind: file | fileT | idx | idxT | mem)
#   environment: BASE, NEW (Esa build directories), ABAPP (the built client), CONT (the directory with work_<input>/out.sesa), PAT (the directory with pat_<input>.txt), ROWS (where the rows go), D (the bench dir, default /tmp/stbench)
set -u
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 DOTNET_gcServer=0
unset DOTNET_GCHeapHardLimit
HERE=$(cd "$(dirname "$0")" && pwd)
D=${D:-/tmp/stbench}; BASE=${BASE:?}; NEW=${NEW:?}; ABAPP=${ABAPP:?}; CONT=${CONT:?}; PAT=${PAT:?}; ROWS=${ROWS:?}
PROCS=${PROCS:-5}; PASSES=${PASSES:-5}; SECS=${SECS:-1.0}
mode=${1:?null|ab}; shift
mkdir -p "$ROWS"
for cell in "$@"; do
  IFS=: read -r input kind th <<< "$cell"
  case $kind in
    mem) suffix="mem:$D/data/$input.txt" ;;
    file|fileT|idx|idxT) suffix="$kind:$CONT/work_$input/out.sesa" ;;
    *) echo "unknown kind $kind" >&2; exit 2 ;;
  esac
  a="$BASE|$suffix"; b="$NEW|$suffix"; [ "$mode" = null ] && b="$BASE|$suffix"
  echo "$(date -u +%H:%M:%SZ) cell $mode $cell load $(cut -d' ' -f1 /proc/loadavg)"
  ABAPP="$ABAPP" "$HERE/abapp/run_ab.sh" "$ROWS/${mode}_${input}_${kind}_t${th}.txt" "$PROCS" "$PAT/pat_$input.txt" "$PASSES" "$th" "$SECS" "$a" "$b" "$mode-$input-$kind-t$th"
done
