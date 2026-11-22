#!/bin/bash
# The in-process A/B client of work package 1b.9g, several processes in a row, rows appended to a file.
#   usage: run_ab.sh <rows file> <processes> <pattern file> <passes> <threads> <seconds> <specA> <specB> [label]
#   spec = <dir>|mem:<text file> | <dir>|file:<container> | <dir>|idx:<container> ("fileT:", "idxT:" = Trusted access); <dir> holds a build of SuffixTree.Esa.dll and its dependencies.
#   environment: ABAPP (the directory of the built client, default ./abapp_bin next to this script).
# Every process is one run of the client (4 warm-up rounds, then <passes> A B B A blocks per operation, every pass <seconds> long with <threads> threads); a line "# <label> proc <i>" precedes its rows.
set -u
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
HERE=$(cd "$(dirname "$0")" && pwd)
ABAPP=${ABAPP:-$HERE/abapp_bin}
out=$1; procs=$2; pat=$3; passes=$4; th=$5; secs=$6; A=$7; B=$8; label=${9:-run}
for p in $(seq 1 "$procs"); do
  echo "# $label proc $p" >> "$out"
  sync
  dotnet "$ABAPP/abapp.dll" "$pat" "$passes" "$th" "$secs" "$A" "$B" | grep -v '^sink' >> "$out"
done
