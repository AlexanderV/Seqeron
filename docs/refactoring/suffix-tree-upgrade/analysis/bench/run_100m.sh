#!/usr/bin/env bash
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $B
export DOTNET_ROOT=$HOME/.dotnet PATH=$HOME/.dotnet:$PATH
rm -f trees/dna_50m.st
python3 lib/benchlib.py measure --raw runs/raw.tsv --data-dir $B --repeats 1 --timeout 600 --tool seqeron-st-persistent \
  --version "Seqeron@242318429a5e net10.0 Release" --input dna_100m --phase build -- ./seqeron/out/SeqeronBench pers-build data/dna_100m.txt trees/dna_100m.st
ls -la trees > runs/trees_after_100m.txt
rm -f trees/dna_100m.st*
echo R100_DONE
