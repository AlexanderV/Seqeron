#!/usr/bin/env bash
# Session 1 of work package 1b.9: the libsais control at the start, the in-memory build rows (A and B of the brief), the in-memory query throughput, the matching operations and the Open of a container.
# One group after the other, nothing else running. Containers of dna_10m and prot_1m must exist (run_1b9_gate.sh container ...). usage: run_1b9_session1.sh   (log on stdout)
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$B/run_1b9_gate.sh; OUT=${OUT:-/tmp/stbench/b9}
stamp() { echo "== $(date -u +%H:%M:%SZ) $* | $(cut -d' ' -f1-3 /proc/loadavg)"; }
stamp "start session 1"; uptime
stamp drift-start;            $G drift start
stamp build-dna_10m;          $G build dna_10m 5
stamp steady-dna_10m;         $G steady dna_10m 3
stamp entry-dna_10m;          $G entry dna_10m 5
for x in dna_1m rep_10m prot_1m dna_20m dna_50m; do stamp build-$x; $G build $x 5; done
for x in dna_10m dna_1m rep_10m prot_1m; do stamp query-$x; $G query $x 5; done
for op in mem mum anchors; do stamp nav-dna_10m-$op; $G nav dna_10m $op 5; done
stamp nav-rep_10m-mem;        $G nav rep_10m mem 5
stamp nav-dna_10m-mem-ukkonen; $G nav dna_10m mem 3 ukkonen
stamp open-dna_10m-warmcache; $G open $OUT/work_dna10/out.sesa dna_10m 5
stamp open-dna_10m-cold;      $G open $OUT/work_dna10/out.sesa dna_10m 5 cold
stamp open-prot_1m-warmcache; $G open $OUT/work_prot1/out.sesa prot_1m 5
stamp open-prot_1m-cold;      $G open $OUT/work_prot1/out.sesa prot_1m 5 cold
stamp drift-mid;              $G drift mid
stamp "end session 1"; uptime
