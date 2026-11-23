#!/usr/bin/env bash
# Session 1b of work package 1b.9 (after session 2, the same quiet machine): the instrument check of the cold build (run_1b9_instr.sh: build / lean / entry in alternating order), the `lean` build rows of every input
# (nothing read or counted before the first stage), the Open of the containers with nothing before the first Open (warm file cache and dropped page cache), then session 2b (abapp: the facade against the index behind it).
# usage: run_1b9_session1b.sh   (log on stdout)
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$B/run_1b9_gate.sh; OUT=${OUT:-/tmp/stbench/b9}
stamp() { echo "== $(date -u +%H:%M:%SZ) $* | $(cut -d' ' -f1-3 /proc/loadavg)"; }
stamp "start session 1b"; uptime
stamp instrument-check;       $B/run_1b9_instr.sh 7
for x in dna_10m dna_1m rep_10m prot_1m dna_20m dna_50m; do stamp lean-$x; $G lean $x 5; done
stamp open-dna_10m-warm;      $G open $OUT/work_dna10/out.sesa dna_10m 5
stamp open-dna_10m-cold;      $G open $OUT/work_dna10/out.sesa dna_10m 5 cold
stamp open-prot_1m-warm;      $G open $OUT/work_prot1/out.sesa prot_1m 5
stamp open-prot_1m-cold;      $G open $OUT/work_prot1/out.sesa prot_1m 5 cold
stamp session-2b;             $B/run_1b9_session2b.sh
stamp "end session 1b"; uptime
