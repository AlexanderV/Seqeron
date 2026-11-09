#!/usr/bin/env bash
# Prototype P2 (work package 1b.4a): end-to-end MEM, 100 kbp query with 1 % substitutions, minLength 20 (the memq100k query of the benchmark generator).
# Fresh process per run, raw rows appended to p2_mem_raw.tsv.   usage: run_p2.sh <bench dir> <input name> <runs> [extra mem options...]
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$B/../../../../.." && pwd)"
export DOTNET_ROOT=${DOTNET_ROOT:-$HOME/.dotnet} PATH=$HOME/.dotnet:$PATH DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_gcServer=0
DLL=$ROOT/apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll
D=$1; X=$2; R=$3; shift 3
RAW=${P2_RAW:-$B/p2_mem_raw.tsv}
COMMIT=$(git -C "$ROOT" rev-parse --short=10 HEAD)$(git -C "$ROOT" diff --quiet HEAD -- src apps 2>/dev/null || echo -dirty)
for r in $(seq $(( ${RUN_BASE:-0} + 1 )) $(( ${RUN_BASE:-0} + R ))); do
  echo "# process input=$X run=$r commit=$COMMIT load=$(cut -d' ' -f1 /proc/loadavg) $(date -u +%FT%TZ) options: $*" >> "$RAW"
  dotnet "$DLL" mem --input "$D/data/$X.txt" --name "$X" --queries "$D" "$@" >> "$RAW" || echo "# process input=$X run=$r EXIT $?" >> "$RAW"
done
