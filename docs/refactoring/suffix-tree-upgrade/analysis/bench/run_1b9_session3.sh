#!/usr/bin/env bash
# Session 3 of work package 1b.9: the file-backed (external-memory) pipeline of S1b on the production chain (analysis/s1b-handoff/s1b.sh): the container build at 10^7 (dna, genome, prot, rep) and 10^8 (dna, and genome / prot / rep when
# EXTRA=1), wall time, RssAnon by the independent sampler, the container size, the SHA-256 against reference_hashes.tsv, the Structure level, the query set, the 11 gates of the checker. dna at 10^7 and 10^8 are built twice more
# (--rebuild) for the spread of the build time. Output: the full chain logs under $OUT/log, the key lines on stdout (results_1b9_chain.txt of the repository is made from them).
# usage: run_1b9_session3.sh   Needs the builds of the Bench app and of EmProto under $OUT/bin (see run_1b9_gate.sh) and about 3 GB of disk at 10^8.
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; OUT=${OUT:-/tmp/stbench/b9}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
unset DOTNET_GCHeapHardLimit
export S1B_BENCH=$OUT/bin/bench/SuffixTree.Esa.Bench.dll S1B_EMPROTO=$OUT/bin/emproto/SuffixTree.EmProto.dll S1B_HEARTBEAT=0
S=$B/../s1b-handoff/s1b.sh; DATA=$OUT/s1bdata; mkdir -p "$DATA"
stamp() { echo "== $(date -u +%H:%M:%SZ) $* | $(cut -d' ' -f1-3 /proc/loadavg)"; }
keys() { grep -E "RESULT|SUMMARY|PASS|FAIL|WARN|total_s|ph1_s|ph23_s|ph45_s|bwt_s|psi_s|file_size|scratch_peak|rssanon_peak|RssAnon|sha256|wall|done in|equals" | cut -c1-210; }
chain() { # <kind> <n>
  stamp "chain $1 $2"; $S chain --n "$2" --kind "$1" --data "$DATA" > "$OUT/log/chain_$1_$2.log" 2>&1; echo "exit $?"; keys < "$OUT/log/chain_$1_$2.log"; }
repeat() { # <kind> <n> <k>: build again (--rebuild) and hash, the lines of the build only
  for k in $(seq 1 "$3"); do stamp "rebuild $1 $2 #$k"; $S build --n "$2" --kind "$1" --data "$DATA" --rebuild > "$OUT/log/rebuild_$1_$2_$k.log" 2>&1; echo "exit $?"
    keys < "$OUT/log/rebuild_$1_$2_$k.log"; $S hash --n "$2" --kind "$1" --data "$DATA" 2>&1 | cut -c1-120; done; }
stamp "start session 3"; uptime
[ "${SKIP_OTHERS:-0}" = 1 ] && { chain dna 10m; repeat dna 10m 2; $S clean --n 10m --kind dna --data "$DATA" > /dev/null 2>&1; chain dna 100m; repeat dna 100m 2; $S clean --n 100m --kind dna --data "$DATA" > /dev/null 2>&1; stamp "end session 3 (dna only)"; exit 0; }
chain dna 10m; repeat dna 10m 2; $S clean --n 10m --kind dna --data "$DATA" > /dev/null 2>&1
for k in genome prot rep; do chain $k 10m; $S clean --n 10m --kind $k --data "$DATA" > /dev/null 2>&1; done
chain dna 100m; repeat dna 100m 2; $S clean --n 100m --kind dna --data "$DATA" > /dev/null 2>&1
[ "${ONLY_DNA:-0}" = 1 ] && { stamp "end session 3"; exit 0; }
if [ "${EXTRA:-0}" = 1 ]; then for k in genome prot rep; do chain $k 100m; $S clean --n 100m --kind $k --data "$DATA" > /dev/null 2>&1; done; fi
stamp "end session 3"; uptime; df -h /tmp | tail -1
