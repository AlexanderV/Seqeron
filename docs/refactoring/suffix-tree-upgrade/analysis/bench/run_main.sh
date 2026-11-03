#!/usr/bin/env bash
# Drives the Seqeron (C#) and Go measurements.  usage: run_main.sh <group> <input> [repeats]
#   groups: go | seqmem | seqpers
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DOTNET_ROOT=$HOME/.dotnet PATH=$HOME/.dotnet:/usr/local/go/bin:$PATH
BL="python3 $B/lib/benchlib.py"
RAW=$B/runs/raw.tsv
SB=$B/seqeron/out/SeqeronBench
GB=$B/go/gobench
SV="Seqeron@242318429a5e net10.0 Release"
GV="go1.24.7 index/suffixarray"
group=$1; X=$2; R=${3:-3}
cd $B
TMO=${TMO:-0}
m() { $BL measure --raw $RAW --data-dir $B --repeats $R --timeout $TMO "$@"; }
case $group in
  go)
    m --tool go-suffixarray --version "$GV" --input $X --phase all -- $GB all data/$X.txt queries/$X
    m --tool go-suffixarray --version "$GV" --input $X --phase serialize --op-prefix ser_ -- $GB write data/$X.txt trees/$X.go.idx
    m --tool go-suffixarray --version "$GV" --input $X --phase load_cold --op-prefix cold_ --drop-caches -- $GB load trees/$X.go.idx queries/$X
    m --tool go-suffixarray --version "$GV" --input $X --phase load_warm --op-prefix warm_ -- $GB load trees/$X.go.idx queries/$X
    rm -f trees/$X.go.idx ;;
  seqmem)
    m --tool seqeron-st-inmem --version "$SV" --input $X --phase all -- $SB mem data/$X.txt queries/$X ;;
  seqpers)
    m --tool seqeron-st-persistent --version "$SV" --input $X --phase build -- $SB pers-build data/$X.txt trees/$X.st
    m --tool seqeron-st-persistent --version "$SV" --input $X --phase query -- $SB pers-query trees/$X.st queries/$X analytics
    m --tool seqeron-st-persistent --version "$SV" --input $X --phase load_cold --op-prefix cold_ --drop-caches -- $SB pers-load trees/$X.st queries/$X
    m --tool seqeron-st-persistent --version "$SV" --input $X --phase load_warm --op-prefix warm_ -- $SB pers-load trees/$X.st queries/$X ;;
esac
