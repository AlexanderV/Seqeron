#!/usr/bin/env bash
# Fresh-mutant audit of work package 1b.3 (ST-S1): runs esa_s1b_analysis_mutants.json against EsaIndex.Analysis.cs and the tests of Esa.Tests/Analysis (default categories only).
# The Release configuration only (the new code uses no Buf<T>). Hangs cannot occur (no child-table walk), but keep a watchdog like run_s1b_audit.sh if mutants are added.
# usage: run_s1b_analysis_audit.sh [first-index last-index]   with ESA_ROOT=<git worktree> to run several slices in parallel: each slice writes results_s1b_analysis_<first>_<last>.tsv
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export MUTANT_RELEASE_ONLY=1
FIRST="${1:-0}"; LAST="${2:-9999}"
python3 - "$FIRST" "$LAST" <<'PY'
import json, sys
a, b = int(sys.argv[1]), int(sys.argv[2])
json.dump(json.load(open('esa_s1b_analysis_mutants.json'))[a:b], open(f'/tmp/esa_s1b_analysis_slice_{a}_{b}.json', 'w'))
PY
python3 run_kernel_mutants.py src/SuffixTree/Algorithms/SuffixTree.Esa/EsaIndex.Analysis.cs "/tmp/esa_s1b_analysis_slice_${FIRST}_${LAST}.json" "results_s1b_analysis_${FIRST}_${LAST}.tsv" "FullyQualifiedName~Esa.Tests.Analysis"
