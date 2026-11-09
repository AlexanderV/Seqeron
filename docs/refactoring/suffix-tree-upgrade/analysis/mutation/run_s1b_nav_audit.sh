#!/usr/bin/env bash
# Fresh-mutant audit of work package 1b.4 (ST-S1): runs the esa_s1b_nav_<set>_mutants.json lists against the new production files and the tests of Esa.Tests/Navigation (default categories).
# Release first, the survivors again in the checked configuration (-p:EsaBounds=true), as run_kernel_mutants.py does.
# usage: ESA_ROOT=<git worktree> run_s1b_nav_audit.sh <set> [<set> ...]   sets: pyramid pack navigator navigation_index matching pinned node navigation_index_b pinned_b navigation_index_c pinned_c pinned_d navigator_b
# Results: results_s1b_nav_<set>.tsv in this directory (id, description, Release result, checked result). Several worktrees can run different sets in parallel. The source file is restored after each set.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
E=src/SuffixTree/Algorithms/SuffixTree.Esa
declare -A FILE=(
  [pyramid]=$E/LcpPyramid.cs [pack]=$E/EsaNavigationPack.cs [navigator]=$E/EsaNavigator.cs [navigation_index]=$E/EsaIndex.Navigation.cs
  [matching]=$E/EsaIndex.Matching.cs [pinned]=$E/PinnedArray.cs [node]=$E/EsaNode.cs
  [navigation_index_b]=$E/EsaIndex.Navigation.cs [pinned_b]=$E/PinnedArray.cs [navigation_index_c]=$E/EsaIndex.Navigation.cs [pinned_c]=$E/PinnedArray.cs [pinned_d]=$E/PinnedArray.cs [navigator_b]=$E/EsaNavigator.cs   # replacements of the five mutants of the first lists that did not compile (analyzer errors)
)
for s in "$@"; do
  python3 run_kernel_mutants.py "${FILE[$s]}" "esa_s1b_nav_${s}_mutants.json" "results_s1b_nav_${s}.tsv" "FullyQualifiedName~Esa.Tests.Navigation"
done
