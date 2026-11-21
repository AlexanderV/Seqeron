#!/bin/bash
# Work package 1b.9d: the cost of the growth of the container alone (`sesa tail` of the Bench application): the three ways of treating a range of S MiB that is added to a file that is then written through 1 MiB mapped
# windows (what PH4 does with the LCP above one byte per entry and with the overflow list): sparse (the behaviour before 1b.9d), zero (the range is reserved by writing zeros: 1b.9d) and worst (the file is created with the
# growth reserved by fallocate: reserving the worst case up front instead). Cold = sync and drop_caches before the run, warm = the page cache as it is. One process per run, the variants alternate inside a round.
# usage: run_1b9d_tail.sh <rows file> <work directory> [rounds, default 3] [sizes in MiB, default "11 100 300 1000"]
set -u
ROWS=$(readlink -f "$1"); W=$2; ROUNDS=${3:-3}; SIZES=${4:-"11 100 300 1000"}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
DLL=${B9D_DLL:-/tmp/stbench/b9d/bin/SuffixTree.Esa.Bench.dll}
mkdir -p "$W"
for r in $(seq 1 "$ROUNDS"); do
  for mb in $SIZES; do
    for cold in true false; do
      VARS="sparse zero worst"
      [ $((r % 2)) = 0 ] && VARS="worst zero sparse"
      for v in $VARS; do
        sync
        dotnet "$DLL" sesa tail --work "$W" --mb "$mb" --variant "$v" --cold "$cold" --label "r$r" >> "$ROWS" 2>> "$ROWS.err"
      done
    done
  done
done
