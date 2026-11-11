#!/usr/bin/env bash
# Fresh-mutant audit of work package 1b.6 (ST-S1): runs the esa_s1b_facade_<set>_mutants.json lists against the new production files of the facade
# (EsaSuffixTree.cs, EsaSuffixTree.Diagnostics.cs, EsaSuffixTree.Navigation.cs, EsaSuffixTreeFactory.cs) with the tests of Esa.Tests/Facade (hand-written classes only in the first pass) and
# IndexConformance (the facade is one of its engines). Release first, the survivors again in the checked configuration (-p:EsaBounds=true); every build and every test run has a watchdog
# (run_index_mutants.py). Second pass: FACADE_FILTER="FullyQualifiedName~Esa.Tests.Facade|FullyQualifiedName~IndexConformance" adds the derived copies of the existing tests (Facade/Ported).
# usage: ESA_ROOT=<git worktree> run_s1b_facade_audit.sh <set> [<start>:<end>] ; sets: core diag nav factory (round 1), core_b diag_b nav_b factory_b (round 2: compiling forms of the mutants that did not compile, and the survivors of round 1 again with the new tests; round 3: surv_core surv_diag surv_nav surv_factory = the survivors of rounds 1 and 2 against the WHOLE default suite, FACADE_FILTER=FullyQualifiedName~SuffixTree.Esa.Tests)
# Results: results_s1b_facade_<set>[_<start>].tsv in this directory. The source file is restored after each run.
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
E=src/SuffixTree/Algorithms/SuffixTree.Esa
declare -A FILE=([core]=$E/EsaSuffixTree.cs [diag]=$E/EsaSuffixTree.Diagnostics.cs [nav]=$E/EsaSuffixTree.Navigation.cs [factory]=$E/EsaSuffixTreeFactory.cs [core_b]=$E/EsaSuffixTree.cs [diag_b]=$E/EsaSuffixTree.Diagnostics.cs [nav_b]=$E/EsaSuffixTree.Navigation.cs [factory_b]=$E/EsaSuffixTreeFactory.cs [surv_core]=$E/EsaSuffixTree.cs [surv_diag]=$E/EsaSuffixTree.Diagnostics.cs [surv_nav]=$E/EsaSuffixTree.Navigation.cs [surv_factory]=$E/EsaSuffixTreeFactory.cs)
s=$1
slice=${2:-}
suffix=${slice:+_${slice/:/-}}
filter=${FACADE_FILTER:-"(FullyQualifiedName~Esa.Tests.Facade&FullyQualifiedName!~Ported)|FullyQualifiedName~IndexConformance"}
out=${FACADE_OUT_PREFIX:-results_s1b_facade}
python3 run_index_mutants.py "${FILE[$s]}" "esa_s1b_facade_${s}_mutants.json" "${out}_${s}${suffix}.tsv" "$slice" "$filter"
