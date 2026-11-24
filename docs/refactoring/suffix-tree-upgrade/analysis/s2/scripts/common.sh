#!/usr/bin/env bash
# Shared settings of the S2-0 scripts (ADR-005). Source it; nothing here runs a job.
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
unset DOTNET_GCHeapHardLimit
S2="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"           # analysis/s2
ROOT="$(cd "$S2/../../../../.." && pwd)"                         # repository root
OUT=${S2_OUT:-/tmp/stbench/s20}                                  # scratch (binaries, texts, containers); never inside the repository
D=${SEQERON_BENCH_DIR:-/tmp/stbench}                             # the S1 bench inputs (data/, queries/); read only
LIBSAIS_BIN=${LIBSAIS_BIN:-$D/b9/sa_bench_libsais}               # built as in analysis/bench/run_1b9_gate.sh
WB="dotnet $OUT/bin/wb/WidthBench.dll"
mkdir -p "$OUT/bin" "$OUT/log" "$S2/results"
quiet() { # refuse to time anything on a loaded machine
  local l; l=$(cut -d' ' -f1 /proc/loadavg)
  echo "load $(cut -d' ' -f1-3 /proc/loadavg) at $(date -u +%H:%M:%SZ)"
  ps -eo pid,pcpu,comm --sort=-pcpu | awk 'NR>1 && $2>20 {print "busy: "$0}' | head -3
}
