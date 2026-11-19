#!/bin/bash
# Cost of the reader count of work package 1b.8 (outside the shipped code), on the container of dna_10m, warm, in memory and over the file at Checked and Trusted, two instruments:
#   throughput  SesaReaderCostTests (variant throughput of 1b.7f): CountOccurrences of 20-mers, FindAllOccurrences of 8-mers, absent 24-mers, a child walk, one tree per process, seven rounds in alternating tree orders and
#               alternating before / after order, median of 5 repeats inside a process (the protocol of run_1b7f_g3.sh); BEFORE = the head of 1b.7f (no reader count), AFTER = every public member of EsaSuffixTree counts itself
#   in-process  FacadeGuardCostTests: the facade call against the call of the index behind it, in alternating batches of one process, and the member IsEmpty against the read of the index length (the count itself)
#   usage: run_1b8_dispose_cost.sh <rows file> <dir of the BEFORE build of SuffixTree.Esa.Tests> <dir of the AFTER build> [rounds of the in-process instrument, default 5]
set -u
ROWS=$(readlink -f "$1")
BEFORE=$(readlink -f "$2")
AFTER=$(readlink -f "$3")
INPROC=${4:-5}
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 DOTNET_gcConcurrent=0 DOTNET_TieredPGO=1
W=${G3_WORK:-/tmp/stbench/g3f}
TXT=/tmp/stbench/data/dna_10m.txt
FILE=$W/dna_10m.sesa
[ -f "$FILE" ] || { echo "missing $FILE (build it with run_1b7f_g3.sh or the SesaBuildCostTests container variant)"; exit 1; }
[ "${TP_ROUNDS:-7}" -gt 0 ] && printf '#throughput: label tree count20_per_s locate8_per_s count24absent_per_s childwalk_nodes_per_s repeats\n' > "$ROWS"
printf '#guard: label tree count20_facade_ns count20_index_ns count20_facade_minus_index_ns count20_overhead_pct locate8_facade_ns locate8_index_ns locate8_facade_minus_index_ns locate8_overhead_pct isempty_ns indexlength_ns guard_ns count20_4threads_facade_ns count20_4threads_index_ns count20_4threads_diff_ns count20_4threads_overhead_pct\n' >> "$ROWS"
runtp() { # variant dir round tree
  sync
  SESA_ROWS=$ROWS SESA_LABEL=$1.r$3 SESA_VARIANT=throughput SESA_TREES=$4 SESA_FILE=$FILE SESA_INPUT=$TXT SESA_REPEATS=5 \
    dotnet test "$2/SuffixTree.Esa.Tests.dll" --nologo --filter "FullyQualifiedName~SesaReaderCostTests" > "$W/dispose_cost_last.log" 2>&1 || { echo "failed: $1 $4"; tail -20 "$W/dispose_cost_last.log"; }
}
ROUND=0
TP=${TP_ROUNDS:-7}   # TP_ROUNDS=0 skips the throughput instrument (its rows were taken by an earlier run of this script and are appended to the same rows file)
for ORDER in "memory checked trusted" "trusted checked memory" "checked memory trusted" "memory trusted checked" "trusted memory checked" "checked trusted memory" "memory checked trusted"; do
  ROUND=$((ROUND + 1))
  [ "$ROUND" -le "$TP" ] || continue
  for TREE in $ORDER; do
    if [ $((ROUND % 2)) -eq 1 ]; then PAIR="before after"; else PAIR="after before"; fi
    for V in $PAIR; do
      if [ "$V" = before ]; then DIR=$BEFORE; else DIR=$AFTER; fi
      runtp "$V" "$DIR" "$ROUND" "$TREE"
    done
  done
done
for K in $(seq 1 "$INPROC"); do
  if [ $((K % 2)) -eq 1 ]; then PAIR="before after"; else PAIR="after before"; fi
  for V in $PAIR; do
    if [ "$V" = before ]; then DIR=$BEFORE; else DIR=$AFTER; fi
    sync
    SESA_ROWS=$ROWS SESA_LABEL=$V.k$K SESA_FILE=$FILE SESA_INPUT=$TXT \
      dotnet test "$DIR/SuffixTree.Esa.Tests.dll" --nologo --filter "FullyQualifiedName~FacadeGuardCostTests" > "$W/dispose_cost_last.log" 2>&1 || { echo "failed: guard $V"; tail -20 "$W/dispose_cost_last.log"; }
  done
done
echo done
