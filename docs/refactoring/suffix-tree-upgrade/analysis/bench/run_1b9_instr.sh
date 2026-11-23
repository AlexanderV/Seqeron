#!/usr/bin/env bash
# The instrument check of work package 1b.9: the cold build of dna_10m through three entries of the same library code in alternating order, <rounds> rounds (default 7), every call a fresh process:
#   build = gate1b9 build (reads /proc/self/status and the GC counters before the first stage), lean = gate1b9 lean (nothing before the build), entry = gate1b9 entry (EsaSuffixTreeFactory.CreateEsa).
# Raw rows -> runs_esa/raw_1b9_instr.tsv (benchlib format; op prefixes lean_ and entry_ are written by the groups). usage: run_1b9_instr.sh [rounds]
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$B/run_1b9_gate.sh
export RAW=$B/runs_esa/raw_1b9_instr.tsv
for r in $(seq 1 "${1:-7}"); do
  echo "== round $r $(date -u +%H:%M:%SZ) $(cut -d' ' -f1-3 /proc/loadavg)"
  case $((r % 3)) in
    0) order="build lean entry" ;; 1) order="lean entry build" ;; *) order="entry build lean" ;;
  esac
  for g in $order; do $G $g dna_10m 1 2>&1 | tail -1; done
done
