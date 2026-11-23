#!/usr/bin/env bash
# The test sets of work package 1b.9 that had not been run on the code of this head: the clean Release build of the solution, SuffixTree.Tests (the Ukkonen and the core tests; the two explicit Performance
# fixtures alone afterwards, on the idle machine), SuffixTree.Mcp.Core.Tests, the Vectors category of SuffixTree.Esa.Tests (the inputs under SEQERON_BENCH_DIR, the libsais vectors in analysis/sa_vectors.tsv) and
# the 10-minute budget run of the 1b.8 differential harness (analysis/em/run_1b8_fuzz.sh). One step after the other; never while a timing run is going. Full logs under $OUT/log, the summary lines on stdout.
# usage: run_1b9_tests.sh [fuzz seed]   (default 20261123)
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$B/../../../../.." && pwd)"; OUT=${OUT:-/tmp/stbench/b9}
SEED=${1:-20261123}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
unset DOTNET_GCHeapHardLimit
cd "$ROOT"
stamp() { echo "== $(date -u +%H:%M:%SZ) $* | $(cut -d' ' -f1-3 /proc/loadavg)"; }
sum() { grep -aE "Passed!|Failed!|Skipped!|error|Warning\(s\)|Error\(s\)|Time Elapsed|Total tests|Passed:|Failed:|Skipped:" "$1" | tail -8 | cut -c1-200; }
stamp "start tests; head $(git rev-parse --short=12 HEAD); git diff 12495c3e..HEAD -- src apps: $(git diff 12495c3e..HEAD --stat -- src apps | wc -l) lines"
stamp "build: dotnet build -c Release --no-incremental Seqeron.sln"
dotnet build Seqeron.sln -c Release --no-incremental --nologo > "$OUT/log/tests_build.log" 2>&1; echo "exit $?"; sum "$OUT/log/tests_build.log"
T=tests/SuffixTree
stamp "SuffixTree.Tests default (Release)"
dotnet test $T/SuffixTree.Tests/SuffixTree.Tests.csproj -c Release --no-build --nologo > "$OUT/log/tests_ukkonen.log" 2>&1; echo "exit $?"; sum "$OUT/log/tests_ukkonen.log"
stamp "SuffixTree.Tests Category=Performance (the two explicit fixtures, idle machine)"
dotnet test $T/SuffixTree.Tests/SuffixTree.Tests.csproj -c Release --no-build --nologo --filter "Category=Performance" > "$OUT/log/tests_ukkonen_perf.log" 2>&1; echo "exit $?"; sum "$OUT/log/tests_ukkonen_perf.log"
stamp "SuffixTree.Mcp.Core.Tests (Release)"
dotnet test $T/SuffixTree.Mcp.Core.Tests/SuffixTree.Mcp.Core.Tests.csproj -c Release --no-build --nologo > "$OUT/log/tests_mcp.log" 2>&1; echo "exit $?"; sum "$OUT/log/tests_mcp.log"
stamp "SuffixTree.Esa.Tests Category=Vectors (Release; inputs in SEQERON_BENCH_DIR=$SEQERON_BENCH_DIR)"
dotnet test $T/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --no-build --nologo --filter "Category=Vectors" > "$OUT/log/tests_vectors.log" 2>&1; echo "exit $?"; sum "$OUT/log/tests_vectors.log"
stamp "1b.8 differential harness, 600 s, master seed $SEED, 4 workers, FUZZ_HANG_S 900"
FUZZ_HANG_S=900 "$B/../em/run_1b8_fuzz.sh" "$B/../em/rows_1b9_fuzz_10min.tsv" 600 "$SEED" 4
grep -E "^# (lane .* done|test process|Passed|Failed)" "$B/../em/rows_1b9_fuzz_10min.tsv" | cut -c1-160
stamp "end tests"
