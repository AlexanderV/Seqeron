#!/usr/bin/env bash
# Session 2b of work package 1b.9: the gate of 1b.9g.1 as it was defined there (the facade over the container against the index behind it, four threads, Checked), on the benchmark's pattern sets, so that the figure
# of session 2 (the file-backed facade against the IN-MEMORY tree, as the brief of 1b.9 asks) can be read next to it. One order (A = the index behind the facade, B = the facade). Rows -> results_1b9_abapp.txt.
# usage: run_1b9_session2b.sh   (the containers and the pattern files of session 2 must exist)
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$B/run_1b9_gate.sh; OUT=${OUT:-/tmp/stbench/b9}; ESA=${ESA_BIN:-$OUT/bin/esa}; ROWS=$B/results_1b9_abapp.txt
stamp() { echo "== $(date -u +%H:%M:%SZ) $* | $(cut -d' ' -f1-3 /proc/loadavg)"; }
stamp "start session 2b"
for x in dna_10m prot_1m; do
  c=$OUT/work_dna10/out.sesa; [ "$x" = prot_1m ] && c=$OUT/work_prot1/out.sesa
  stamp "$x.4t.idx_file"; $G ab "$ROWS" 3 $OUT/pat_$x.txt 9 4 1 "$ESA|idx:$c" "$ESA|file:$c" "$x.4t.idx_file"
done
stamp "end session 2b"
