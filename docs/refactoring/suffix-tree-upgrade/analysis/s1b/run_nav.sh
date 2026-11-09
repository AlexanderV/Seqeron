#!/usr/bin/env bash
# Work package 1b.4: end-to-end operations of the production navigator (MEM, MUM, anchors of the 100 kbp query with 1 % substitutions at minLength 20; the LCS family; ns per link; the MEM time split).
# Fresh process per run, raw rows appended to nav_raw.tsv.   usage: run_nav.sh <bench dir> <input name> <op> <runs> [extra nav options...]
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$B/../../../../.." && pwd)"
export DOTNET_ROOT=${DOTNET_ROOT:-$HOME/.dotnet} PATH=$HOME/.dotnet:$PATH DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_gcServer=0
DLL=${NAV_DLL:-$ROOT/apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll}
D=$1; X=$2; OP=$3; R=$4; shift 4
RAW=${NAV_RAW:-$B/nav_raw.tsv}
COMMIT=$(git -C "$ROOT" rev-parse --short=10 HEAD)$(git -C "$ROOT" diff --quiet HEAD -- src apps 2>/dev/null || echo -dirty)
for r in $(seq $(( ${RUN_BASE:-0} + 1 )) $(( ${RUN_BASE:-0} + R ))); do
  echo "# process input=$X op=$OP run=$r commit=$COMMIT load=$(cut -d' ' -f1 /proc/loadavg) $(date -u +%FT%TZ) options: $*" >> "$RAW"
  dotnet "$DLL" nav --input "$D/data/$X.txt" --name "$X" --queries "$D" --op "$OP" "$@" >> "$RAW" || echo "# process input=$X op=$OP run=$r EXIT $?" >> "$RAW"
done
