#!/usr/bin/env bash
# Fresh-mutant audit of work packages 1b.0-1b.2 (ST-S1): runs every mutant list of this directory against the file it mutates and the test classes that own it.
# Hangs: a mutant that makes the descent or the scan loop forever never ends; run a watchdog next to this script that kills a test host older than 180 s (counted as killed).
# usage: run_s1b_audit.sh [set-name ...]    (default: all sets). Results: results_s1b_<set>.tsv (id, description, Release result). The source file is restored after each set.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
ROOT="$(cd ../../../../.. && pwd)"
export MUTANT_RELEASE_ONLY=1
E=src/SuffixTree/Algorithms/SuffixTree.Esa
declare -A FILE=(
  [lcp_interval_scan]=$E/LcpIntervalScan.cs [child_layouts]=$E/ChildLayouts.cs [child_layouts_tagged]=$E/ChildLayouts.cs [child_table]=$E/ChildTable.cs
  [leftmost_occurrence]=$E/LeftmostOccurrence.cs [traverse]=$E/EsaIndex.Traverse.cs [search]=$E/EsaIndex.Search.cs
  [contract_locus]=src/SuffixTree/Algorithms/SuffixTree.Core/Indexing/IndexLocus.cs
  [contract_adapter]=tests/SuffixTree/SuffixTree.Esa.Tests/Contract/UkkonenIndexAdapter.cs
)
declare -A FILTER=(
  [lcp_interval_scan]="FullyQualifiedName~Esa.Tests.Tree" [child_layouts]="FullyQualifiedName~Esa.Tests.Tree" [child_layouts_tagged]="FullyQualifiedName~Esa.Tests.Tree" [child_table]="FullyQualifiedName~Esa.Tests.Tree"
  [leftmost_occurrence]="FullyQualifiedName~Esa.Tests.Tree" [traverse]="FullyQualifiedName~Esa.Tests.Tree" [search]="FullyQualifiedName~Esa.Tests.Search.EsaSearchTests"
  [contract_locus]="FullyQualifiedName~Esa.Tests.Contract" [contract_adapter]="FullyQualifiedName~Esa.Tests.Contract"
)
sets=("$@"); [ ${#sets[@]} -eq 0 ] && sets=(lcp_interval_scan child_layouts child_layouts_tagged child_table leftmost_occurrence traverse search contract_locus contract_adapter)
for s in "${sets[@]}"; do
  python3 run_kernel_mutants.py "${FILE[$s]}" "esa_s1b_${s}_mutants.json" "results_s1b_${s}.tsv" "${FILTER[$s]}"
done
