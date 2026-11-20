#!/bin/bash
# The budgeted run of the full differential fuzz of work package 1b.8 (gate G2), outside the shipped code: the explicit class FullDifferentialFuzzTests of SuffixTree.Esa.Tests in the Release build, the whole matrix
# (the generator families x alphabets x size bands over the in-memory and the file-backed engines against the definitions, the Ukkonen tree and the v6 tree; the container corruption fuzz; the shared-index test), one master
# seed, a time budget (default 1800 s) split among the lanes by share (FullFuzzOptions.DefaultLanes). The machine must be otherwise idle. The report (lanes, per-engine counts and throughput, the corruption and concurrency
# counts, the calls of every member) and the progress log are folded into the rows file: lines of the report as they are, progress lines behind a '#'.
#   usage: run_1b8_fuzz.sh <rows file> [budget seconds, default 1800] [master seed, default 1] [threads, default the number of cores]   (environment FUZZ_HANG_S: the watchdog limit of a case in seconds, default 300; 1b.9a: 900 for the xl band)
set -u
ROWS=$(readlink -f "$1")
BUDGET=${2:-1800}
SEED=${3:-1}
THREADS=${4:-$(nproc)}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
PROJ=$ROOT/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj
W=${FUZZ_WORK:-/tmp/fuzz1b8}
mkdir -p "$W"
rm -f "$W/report.txt" "$W/report.txt.log"
rm -rf "$W/journal"
{
  printf '# run_1b8_fuzz.sh: head %s, budget %s s, master seed %s, threads %s, started %s\n' "$(git -C "$ROOT" rev-parse HEAD)" "$BUDGET" "$SEED" "$THREADS" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '# uptime before: %s\n' "$(uptime)"
} > "$ROWS"
START=$(date +%s)
ESA_FUZZ_SEED=$SEED ESA_FUZZ_BUDGET_S=$BUDGET ESA_FUZZ_THREADS=$THREADS ESA_FUZZ_HANG_S=${FUZZ_HANG_S:-300} ESA_FUZZ_REPORT=$W/report.txt ESA_FUZZ_JOURNAL=$W/journal \
  dotnet test "$PROJ" -c Release --no-build --nologo --filter "FullyQualifiedName~FullDifferentialFuzzTests" > "$W/test.out" 2>&1
CODE=$?
END=$(date +%s)
{
  printf '# test process exit code %s, wall clock %s s, ended %s\n' "$CODE" "$((END - START))" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '# uptime after: %s\n' "$(uptime)"
  grep -E "Passed!|Failed!|HANG|crashed" "$W/test.out" | sed 's/^/# /'
  [ -f "$W/report.txt.log" ] && sed 's/^/# /' "$W/report.txt.log"
  [ -f "$W/report.txt" ] && cat "$W/report.txt"
} >> "$ROWS"
echo "exit $CODE, rows in $ROWS"
