#!/usr/bin/env bash
# Prototype P1 (work package 1b.4a): suffix-link kernels. Fresh process per run, raw `R<TAB>...` rows appended to p1_links_raw.tsv with a `# process` header line.
#   usage: run_p1.sh <bench dir> <input name> <runs> [extra links options...]
#   inputs: dna_10m | rep_10m | prot_1m (files <bench dir>/data/<name>.txt)  or  a^n | fibonacci | thue-morse (--synthetic <name>:1000000)
# Quiet machine only (uptime, ps). The kernel order is rotated by run so that no kernel always runs first.
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$B/../../../../.." && pwd)"
export DOTNET_ROOT=${DOTNET_ROOT:-$HOME/.dotnet} PATH=$HOME/.dotnet:$PATH DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_gcServer=0
DLL=$ROOT/apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll
D=$1; X=$2; R=$3; shift 3
RAW=${P1_RAW:-$B/p1_links_raw.tsv}
COMMIT=$(git -C "$ROOT" rev-parse --short=10 HEAD)$(git -C "$ROOT" diff --quiet HEAD -- src apps 2>/dev/null || echo -dirty)
case "$X" in
  dna_10m|rep_10m|prot_1m|dna_1m|dna_20m|dna_50m) SRC=(--input "$D/data/$X.txt" --name "$X") ;;
  *) SRC=(--synthetic "${X}") ;;
esac
for r in $(seq 1 "$R"); do
  echo "# process input=$X run=$r commit=$COMMIT load=$(cut -d' ' -f1 /proc/loadavg) $(date -u +%FT%TZ) options: $*" >> "$RAW"
  dotnet "$DLL" links "${SRC[@]}" --rotate $(( (r - 1) * 2 )) "$@" >> "$RAW" || echo "# process input=$X run=$r EXIT $?" >> "$RAW"
done
