#!/bin/bash
# Work package 1b.9d: the in-process A/B of the growth of the container on whole builds (`sesa ab` of the Bench application): the production build, default parameters, four threads, the same process alternating A (the growth is
# sparse: the behaviour before the package) and B (the growth is reserved by writing zeros), <pairs> pairs, the order alternating from pair to pair; cold = sync and drop_caches before every build, warm = the page cache as it is.
# usage: run_1b9d_ab.sh <rows file> <work directory> <pairs> <cold true|false> <label=input file>...
set -u
ROWS=$(readlink -f "$1"); W=$2; PAIRS=$3; COLD=$4; shift 4
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
DLL=${B9D_DLL:-/tmp/stbench/b9d/bin/SuffixTree.Esa.Bench.dll}
mkdir -p "$W"
for spec in "$@"; do
  L=${spec%%=*}; F=${spec#*=}
  sync
  dotnet "$DLL" sesa ab --input "$F" --work "$W" --label "$L.$([ "$COLD" = true ] && echo cold || echo warm)" --pairs "$PAIRS" --cold "$COLD" --threads 4 >> "$ROWS" 2>> "$ROWS.err"
done
