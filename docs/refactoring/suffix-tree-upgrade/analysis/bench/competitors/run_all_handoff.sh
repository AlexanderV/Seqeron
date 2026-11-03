#!/usr/bin/env bash
# Runs every competitor script one after the other (never in parallel: the scripts time multi-threaded throughput).
# usage: BENCH_DIR=/path/to/bench competitors/run_all_handoff.sh        (inputs must exist: see HANDOFF.md section 0)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for s in run_libsais_libdivsufsort.sh run_rustbio.sh run_python.sh run_aligners.sh run_sdsl.sh; do
  echo "=== $s ==="; bash "$HERE/$s" || echo "[$s exited with $?] continuing"
done
echo "all done; assemble with: python3 \$BENCH_DIR/lib/assemble.py \$BENCH_DIR"
