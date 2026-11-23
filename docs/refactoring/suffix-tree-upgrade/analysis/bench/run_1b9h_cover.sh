#!/bin/bash
# Work package 1b.9h, item 4 (the difference cover): the measurements of ADR-004 section 5a "1b.9h", in the order they were made. Nothing else may run on the machine; never build or test while it runs.
# usage: run_1b9h_cover.sh <work dir> <commit before> <commit after>      (the repository checkout is the current directory; paths below are relative to it)
#   1. two builds of SuffixTree.Esa.dll, the commit before and the commit after (git worktrees of the two commits), and the in-process A/B client abcover
#   2. where the time of one cover goes (abcover split), the null experiment (the same build on both sides), the A/B of the two builds (abcover ab), the two covers residue by residue (abcover same)
#   3. the whole build in one process, alternating (the Bench command `sesa resume abhead`, Resume off), on a text of 10^7 symbols, the null experiment first
#   4. the whole build in fresh processes, alternating (the Bench of each commit, `sesa build`), on the same text
set -eu
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
W=$(cd "$1" && pwd); BEFORE=$2; AFTER=$3
REPO=$(git rev-parse --show-toplevel)
for v in before after; do
  c=$BEFORE; [ $v = after ] && c=$AFTER
  [ -d "$W/wt_$v" ] || git worktree add "$W/wt_$v" "$c"
  dotnet build "$W/wt_$v/src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj" -c Release -o "$W/esa_$v" --nologo -v q
  dotnet publish "$W/wt_$v/apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj" -c Release -o "$W/bench_$v" --nologo -v q -p:UseAppHost=false
done
dotnet build "$REPO/docs/refactoring/suffix-tree-upgrade/analysis/bench/abcover/abcover.csproj" -c Release -o "$W/abcover" --nologo -v q
AB="dotnet $W/abcover/abcover.dll"
echo "== where the time goes (the commit before)";            $AB split "$W/esa_before" 1024 20
echo "== null experiment, A = B = the commit before";         $AB ab "$W/esa_before" "$W/esa_before" 1024 20 9
echo "== A/B, A = before, B = after";                         $AB ab "$W/esa_before" "$W/esa_after" 1024 20 9
echo "== the two covers";                                     $AB same "$W/esa_before" "$W/esa_after" 1024
