#!/bin/bash
# work package 1b.9h, item 4 (G4): the containers of the 10^6 and 10^7 cells of reference_hashes.tsv built by the new binaries (g4 part), and the 10^8 DNA chain (the second part); the paths are those of the scratch directory of that package
# G4: the containers of the hand-off reference set at 10^6 and 10^7 (dna, genome, prot, rep) built by the new binaries; the checker compares the SHA-256 with reference_hashes.tsv
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 SEQERON_BENCH_DIR=/tmp/stbench
V=${1:-new}; D=/tmp/stbench/b9h/s1bdata_$V
export S1B_BENCH=/tmp/stbench/b9h/bin_$V/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.dll S1B_EMPROTO=/tmp/stbench/b9h/bin_$V/SuffixTree.EmProto/SuffixTree.EmProto.dll S1B_HEARTBEAT=0
S=/home/user/Seqeron/docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff/s1b.sh
mkdir -p $D
for n in 1m 10m; do for k in dna genome prot rep; do
  echo "== $V $k $n"; $S chain --n $n --kind $k --structure no --queries no --data $D 2>&1 | grep -E "RESULT|equals the reference|container sha256|sha256|differs|FAIL" | head -6 | cut -c1-200
  $S clean --n $n --kind $k --data $D > /dev/null 2>&1
done; done

# the 10^8 DNA chain: one production build, RssAnon against 320 MB, the container hash against the reference
V=${1:-new}; D=/tmp/stbench/b9h/s1bdata_$V; $S chain --n 100m --kind dna --structure no --queries no --data $D
