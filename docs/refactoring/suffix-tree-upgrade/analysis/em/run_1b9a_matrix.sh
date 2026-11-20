#!/bin/bash
# Work package 1b.9a (outside the shipped code): the read-path gap of the file-backed tree. One tree per process, the variants in a ROTATED order from round to round
# (round r starts with the (r-1) mod k-th variant), warm cache, the harness of analysis/em/readbench. Rows go to a TSV file; summarize_1b9a.py reports the medians.
#   usage: run_1b9a_matrix.sh <rows file> <rounds> "<variant list>" <bin dir of ReadBench> <container> <text> [alphabet] [passes]
#   a variant may carry an environment prefix: name=ENV1:VAL1,ENV2:VAL2 (the label is the name); the bin dir of a variant can be overridden with name@bindir
set -u
ROWS=$(readlink -f -m "$1"); ROUNDS=$2; VARIANTS=($3); BIN=$4; FILE=$5; TEXT=$6; ALPHA=${7:-ACGT}; PASSES=${8:-7}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
K=${#VARIANTS[@]}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
printf '#head\t%s\n#row: label variant count20_per_s locate8_per_s absent24_per_s open_s sink smaps...\n' "$(git -C "$ROOT" rev-parse HEAD)" >> "$ROWS"
for ((r = 0; r < ROUNDS; r++)); do
  for ((j = 0; j < K; j++)); do
    spec=${VARIANTS[$(((r + j) % K))]}
    envs=""; bin=$BIN; tag=""
    if [[ "$spec" == *"@"* ]]; then bin=${spec#*@}; spec=${spec%%@*}; tag="@$(basename "$bin")"; fi
    if [[ "$spec" == *"="* ]]; then envs=${spec#*=}; spec=${spec%%=*}; fi
    name=$spec
    envarr=(); IFS=',' read -ra pairs <<< "$envs"; for p in "${pairs[@]}"; do [ -n "$p" ] && envarr+=("${p%%:*}=${p#*:}"); done
    sync
    env "${envarr[@]}" RB_TAG="$tag" dotnet "$bin/SuffixTree.Esa.Bench.dll" "$name" "$FILE" "$TEXT" "$PASSES" "r$((r + 1))" "$ALPHA" >> "$ROWS" || echo "failed: $name round $((r + 1))" >&2
  done
done
