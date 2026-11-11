#!/usr/bin/env bash
# Fresh-mutant audit of work package 1b.5 (ST-S1): runs the esa_s1b_index_<set>_mutants.json lists against the new production files (EsaSuffixIndex*.cs, EsaIndexText.cs) and the tests of Esa.Tests/Contract
# (default categories). Release first, the survivors again in the checked configuration (-p:EsaBounds=true); every build and test run has a watchdog (run_index_mutants.py).
# usage: ESA_ROOT=<git worktree> run_s1b_index_audit.sh <set> [<start>:<end>] ; sets: core tree arrays text core_b tree_b text_b (second round)
# Results: results_s1b_index_<set>[_<start>].tsv in this directory. Several worktrees can run different sets or slices in parallel. The source file is restored after each run.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
E=src/SuffixTree/Algorithms/SuffixTree.Esa
declare -A FILE=([core]=$E/EsaSuffixIndex.cs [tree]=$E/EsaSuffixIndex.Tree.cs [arrays]=$E/EsaSuffixIndex.Arrays.cs [text]=$E/EsaIndexText.cs [core_b]=$E/EsaSuffixIndex.cs [tree_b]=$E/EsaSuffixIndex.Tree.cs [text_b]=$E/EsaIndexText.cs)
s=$1
slice=${2:-}
suffix=${slice:+_${slice/:/-}}
python3 run_index_mutants.py "${FILE[$s]}" "esa_s1b_index_${s}_mutants.json" "results_s1b_index_${s}${suffix}.tsv" "$slice"
