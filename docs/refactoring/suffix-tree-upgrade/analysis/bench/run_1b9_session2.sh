#!/usr/bin/env bash
# Session 2 of work package 1b.9: the in-process A/B client (abapp) on the benchmark's own large pattern sets (100,000 present 20-mers, 1,000 8-mers for locate (4-mers on protein), 100,000 absent 20-mers):
# absolute throughput of the in-memory tree and the file-backed tree (Checked, Trusted) at one and at four threads, the null experiment first, both orders of A and B.
# The in-memory tree is the facade over the built index (EsaSuffixTreeFactory.Create), the file-backed tree the facade over the container (SuffixIndexFactory.Open, Header level).
# usage: run_1b9_session2.sh   (rows -> results_1b9_abapp.txt; the log on stdout). The containers of dna_10m and prot_1m must exist (run_1b9_gate.sh container ...).
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; G=$B/run_1b9_gate.sh; OUT=${OUT:-/tmp/stbench/b9}; D=${D:-/tmp/stbench}; ESA=${ESA_BIN:-$OUT/bin/esa}
ROWS=$B/results_1b9_abapp.txt
stamp() { echo "== $(date -u +%H:%M:%SZ) $* | $(cut -d' ' -f1-3 /proc/loadavg)"; }
stamp "start session 2"; uptime
$G patterns dna_10m $OUT/pat_dna_10m.txt; $G patterns prot_1m $OUT/pat_prot_1m.txt
specs() { # <input> <container>: sets MEM FILE FILET
  MEM="$ESA|mem:$D/data/$1.txt"; FILE="$ESA|file:$2"; FILET="$ESA|fileT:$2"; }
run() { # <label> <threads> <passes> <pattern file> <A> <B>
  stamp "$1"; $G ab "$ROWS" 3 "$4" "$3" "$2" 1 "$5" "$6" "$1"; }
specs dna_10m $OUT/work_dna10/out.sesa; P=$OUT/pat_dna_10m.txt
run dna_10m.4t.null_mem_mem        4 9 $P "$MEM"  "$MEM"
run dna_10m.4t.mem_file            4 9 $P "$MEM"  "$FILE"
run dna_10m.4t.file_mem            4 9 $P "$FILE" "$MEM"
run dna_10m.4t.file_filetrusted    4 9 $P "$FILE" "$FILET"
specs prot_1m $OUT/work_prot1/out.sesa; P=$OUT/pat_prot_1m.txt
run prot_1m.4t.null_mem_mem        4 9 $P "$MEM"  "$MEM"
run prot_1m.4t.mem_file            4 9 $P "$MEM"  "$FILE"
run prot_1m.4t.file_mem            4 9 $P "$FILE" "$MEM"
specs dna_10m $OUT/work_dna10/out.sesa; P=$OUT/pat_dna_10m.txt
run dna_10m.1t.null_mem_mem        1 7 $P "$MEM"  "$MEM"
run dna_10m.1t.mem_file            1 7 $P "$MEM"  "$FILE"
run dna_10m.1t.file_filetrusted    1 7 $P "$FILE" "$FILET"
specs prot_1m $OUT/work_prot1/out.sesa; P=$OUT/pat_prot_1m.txt
run prot_1m.1t.mem_file            1 7 $P "$MEM"  "$FILE"
stamp "end session 2"; uptime
