# own nav A/B of the verification of S2.3 (scratch paths: W=/tmp/stbench/orch33 holds bin_base, bin_base2 (a copy of bin_base), bin_new = builds of apps/SuffixTree.Esa.Bench at 56456a72 and cc047d79)
#!/bin/bash
# own A/B of the in-memory navigator: BASE (56456a72 = S1 head production code) vs NEW (HEAD), plus a null arm (a second copy of BASE); order rotated every round
# usage: nav_ab.sh <op> <rounds> <out.tsv>
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 DOTNET_gcServer=0
unset DOTNET_GCHeapHardLimit
D=/tmp/stbench; W=/tmp/stbench/orch33
op=$1; rounds=$2; out=$3
arms=(base base2 new)
perms=("0 1 2" "1 2 0" "2 0 1" "2 1 0" "1 0 2" "0 2 1")
: > $out
for r in $(seq 1 $rounds); do
  for a in ${perms[$(( (r-1) % 6 ))]}; do
    arm=${arms[$a]}
    txt=$(dotnet $W/bin_$arm/SuffixTree.Esa.Bench.dll nav --input $D/data/dna_10m.txt --name dna_10m --queries $D --op $op --min-length 20 --calls ${CALLS:-10} --ukkonen false 2>&1)
    v=$(echo "$txt" | grep -E "${op}_(warm|steady)_median_s" | tr "\n" "|")
    echo -e "$r\t$arm\t$v" >> $out
  done
done
echo done >> $out
