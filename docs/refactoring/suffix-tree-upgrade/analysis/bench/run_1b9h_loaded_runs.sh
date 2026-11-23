#!/bin/bash
# usage: loaded_runs.sh <runs> <label> <load procs> <filter> ; one dotnet test run per iteration (Release, --no-build) with N busy python processes beside it; writes one line per run and a summary: passed / failed / skipped (an inconclusive test is a skipped one)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
N=$1; LABEL=$2; LOAD=$3; FILT=$4
OUT=/tmp/stbench/b9h/loaded_$LABEL.log; : > $OUT
pids=""; for i in $(seq 1 $LOAD); do python3 -c "
while True: pass" & pids="$pids $!"; done
pass=0; fail=0; skipped=0
for i in $(seq 1 $N); do
  r=$(dotnet test /home/user/Seqeron/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --no-build --filter "$FILT" 2>&1)
  line=$(echo "$r" | grep -E "^(Passed|Failed)!" | head -1 | cut -c1-120)
  sk=$(echo "$line" | sed -n 's/.*Skipped: *\([0-9]*\).*/\1/p'); [ -z "$sk" ] && sk=0
  skipped=$((skipped+sk))
  if echo "$line" | grep -q "^Passed!"; then pass=$((pass+1)); else fail=$((fail+1)); echo "$r" | grep -E "Failed |Error Message" | head -3 | cut -c1-200 >> $OUT; fi
  echo "run $i: $line" >> $OUT
done
kill $pids 2>/dev/null
echo "label $LABEL load $LOAD runs $N passed $pass failed $fail skipped_tests_total $skipped filter $FILT" >> $OUT
