#!/usr/bin/env bash
# Drives the ESA backend measurements (ST-S1 work package 1a.8) with the shared driver lib/benchlib.py: every repeat is a FRESH process.
#   usage: run_esa.sh <group> <data-dir> <input> [repeats]       groups: build | steady | profile
#          run_esa.sh aggregate                                  raw rows -> bench_esa.tsv (median over repeats, min/max in the notes)
#   <data-dir> holds data/<input>.txt made by gen/gen_inputs.py (outside the repository, hashes in MANIFEST.sha256).
#   Build the harness first:  dotnet build apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj -c Release
#   Timing runs only on a quiet machine (uptime, ps); bench.tsv is never touched: raw rows go to runs_esa/raw_esa.tsv, the result to bench_esa.tsv.
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$B/../../../../.." && pwd)"
export DOTNET_ROOT=${DOTNET_ROOT:-$HOME/.dotnet} PATH=$HOME/.dotnet:$PATH DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_gcServer=0
DLL=$ROOT/apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll
BL="python3 $B/lib/benchlib.py"
RAW=$B/runs_esa/raw_esa.tsv
VEC=$B/../sa_vectors.tsv
VARIANT=${ESA_VARIANT:-SaisReference+LcpReference}
SV="Seqeron@$(git -C "$ROOT" rev-parse --short=12 HEAD)$(git -C "$ROOT" diff --quiet HEAD -- src apps 2>/dev/null || echo -dirty) net10.0 Release $VARIANT"
mkdir -p "$B/runs_esa"
group=$1
if [ "$group" = aggregate ]; then
  exec $BL aggregate --raw "$RAW" --out "$B/bench_esa.tsv"
fi
D=$2; X=$3; R=${4:-5}
m() { $BL measure --raw "$RAW" --data-dir "$B" --tool seqeron-esa --version "$SV" --input "$X" "$@"; }
case $group in
  build)   m --repeats "$R" --phase build   -- dotnet "$DLL" build   --input "$D/data/$X.txt" --name "$X" --verify "$VEC" ;;
  steady)  m --repeats "$R" --phase steady  -- dotnet "$DLL" steady  --input "$D/data/$X.txt" --name "$X" --builds 5 ;;
  profile) m --repeats "$R" --phase profile -- dotnet "$DLL" profile --input "$D/data/$X.txt" --name "$X" --builds 3 ;;
  *) echo "unknown group $group" >&2; exit 2 ;;
esac
