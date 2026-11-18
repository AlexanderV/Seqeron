#!/bin/bash
# G5 of work package 1b.7c (outside the shipped code): the cost of reading the child table through the Packed codec against Tagged, on the production files of a 100 Mbp input. One process = ChildCodecReadCostTests.RunOnce
# (loads the text, the suffix array, the LCP files and the production .cld; builds the Tagged table in memory from the same arrays; checks that the production file equals the in-RAM packed table; then four passes of the same
# query mix on both: the first pass of a process is "cold" (JIT, first touch), the others "warm"); the page cache is dropped before every process; the variant that goes first alternates between processes.
#   usage: run_1b7c_g5.sh <rows file> [processes, default 3]
# environment: G5_NAME (default dna_100m), G5_DIR (default /tmp/stbench/c7: <name>.sa, the production outputs <name>.lcp .ovf .cld and the text), G5_TEXT (default /tmp/stbench/e6/dna_100m.txt)
set -u
ROWS=$(readlink -f "$1")
RUNS=${2:-3}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj
NAME=${G5_NAME:-dna_100m}
D=${G5_DIR:-/tmp/stbench/c7}
TEXT=${G5_TEXT:-/tmp/stbench/e6/dna_100m.txt}
printf '#label\tcodec\tphase\tpass\twalk_s\twalk_child_steps\tcount_present_s\tcount_absent_s\tlocate12_s\n' > "$ROWS"
printf '#head\t%s\n' "$(git -C "$ROOT" rev-parse HEAD)" >> "$ROWS"
for K in $(seq 1 "$RUNS"); do
  if [ $((K % 2)) = 1 ]; then FIRST=packed; else FIRST=tagged; fi
  sync; echo 3 > /proc/sys/vm/drop_caches
  G5_FIRST=$FIRST G5_TEXT=$TEXT G5_SA=$D/$NAME.sa G5_LCP=$D/$NAME G5_CLD=$D/$NAME.cld G5_ROWS=$ROWS G5_LABEL="$NAME.proc$K.first_$FIRST" \
    dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~ChildCodecReadCostTests" > "${G5_WORK:-/tmp/stbench/g3c}/g5_last.log" 2>&1 || { echo "run failed: $K"; tail -20 "${G5_WORK:-/tmp/stbench/g3c}/g5_last.log"; }
done
echo g5-done
