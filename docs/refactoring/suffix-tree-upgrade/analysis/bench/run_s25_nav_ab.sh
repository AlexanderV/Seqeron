#!/bin/bash
# the navigator rows of work package S2.5 (the rule of the orchestrator's verification of S2.3: arms base, base2 = a second copy of base (the null) and new, the order rotated every round, one process per arm and round, at least 60 rounds;
# the figure is the per-process minimum of the warm calls and the paired median of the processes; summarize_s23_orch_nav.py reads the rows).
# BASE builds (apps/SuffixTree.Esa.Bench at the S1 head 56456a72 with the one --container patch of this package), NEW (this package): $W/bin_base_c, $W/bin_base2_c (a copy of bin_base_c), $W/bin_new.
# usage: run_s25_nav_ab.sh <op: anchors|mem|mum> <rounds> <out.tsv> <backing: mem|container> <input name: dna_10m|rep_10m|...> [container path]   (CALLS = calls per process, default 40)
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 DOTNET_gcServer=0
unset DOTNET_GCHeapHardLimit
D=${SEQERON_BENCH_DIR:-/tmp/stbench}; W=${S25_DIR:-/tmp/stbench/s25}
op=$1; rounds=$2; out=$3; backing=$4; name=$5; cont=$6
arms=(base base2 new)
dirs=(bin_base_c bin_base2_c bin_new)
perms=("0 1 2" "1 2 0" "2 0 1" "2 1 0" "1 0 2" "0 2 1")
extra=""
if [ "$backing" = "container" ]; then extra="--container $cont"; fi
: > $out
for r in $(seq 1 $rounds); do
  for a in ${perms[$(( (r-1) % 6 ))]}; do
    txt=$(dotnet $W/${dirs[$a]}/SuffixTree.Esa.Bench.dll nav --input $D/data/$name.txt $extra --name $name --queries $D --op $op --min-length 20 --calls ${CALLS:-40} --ukkonen false 2>&1)
    v=$(echo "$txt" | grep -E "${op}_(warm|steady)_median_s|${op}_results|${op}_checksum" | tr "\n" "|")
    echo -e "$r\t${arms[$a]}\t$v" >> $out
  done
done
echo done >> $out
