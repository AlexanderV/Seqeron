# Suffix-tree upgrade: competitor analysis and plan (2026-11-03)

Owner's goal: make the Seqeron suffix tree better, on average, than every competitor. Method: take all competitors,
take all their capabilities, score each intersection in a two-dimensional table, then derive a plan whose projected
average beats all of them. **Status: proposal, not approved. No stage has been started.**

| file | content |
|---|---|
| [COMPETITOR_MATRIX.md](COMPETITOR_MATRIX.md) | 12 competitors + Seqeron: 24 properties (0-4, weighted) and 81 functions (0-2); ranking; gap list |
| [PLAN.md](PLAN.md) | stages S0-S8: scope, oracle, numeric acceptance gates, projected cells, cost, risks, order, hand-off items, decisions |
| [BENCHMARK.md](BENCHMARK.md) | measurements on one machine: Seqeron (in-memory, persistent), Go, libsais, libdivsufsort, sdsl-lite, rust-bio, pydivsufsort, Python suffix-trees, BWA, Bowtie2 |
| `analysis/` | rubric (`SPEC.md`, `FUNCTIONS.md`), per-cell evidence TSVs, Seqeron audit, benchmark harness and raw rows, hand-off package, scripts |

## Headline

| | Seqeron now | best competitor | after the recommended core (S0-S6, S8) | after all stages |
|---|---|---|---|---|
| property mean (24, weighted, 0-4) | 1.80 | 2.48 (ropebwt3) | 2.80 | 3.02 |
| function coverage (81) | 55.6% (39.5% on the tree layer only) | 52.5% (sdsl) | 93.8% | 95.1% |
| first place over 2000 random weightings | 0% | - | 99.5% | 100% |

Measured at 10 Mbp random DNA, one thread: Seqeron builds in 14.6 s (libsais 0.33 s, Go 0.73 s), uses 169 B/char in memory
(56 B/char as persistent file; suffix array 4 B/char), counts 590k 20-mers/s (best 747k), locates 8-mers 23x slower than the best.
The persistent tree finished 50 Mbp (389 s) but not 100 Mbp in 600 s. Reads scale 4.4x on 4 threads.
The projection is conditional on each stage passing its gate; targets are not promises. See PLAN.md section 5 for what it does not achieve.

## Reproduce

```
cd docs/refactoring/suffix-tree-upgrade/analysis
python3 -I merge.py           # property matrix, weighted ranking
python3 -I merge_func.py      # function coverage, group coverage, gap list
python3 plan_model.py         # projection of the stages (writes nothing)
python3 -I build_docs.py      # regenerates ../COMPETITOR_MATRIX.md
```

Benchmark: `analysis/bench/HANDOFF.md` (inputs are regenerated deterministically; `MANIFEST.sha256` lists the hashes measured).
Measured third-party tools were built from open-source repositories only; versions, commit SHAs and licences are in `analysis/bench/SOURCES.txt`.
Versions measured can differ from versions scored in the matrix (for example sdsl-lite 2.1.1 measured, xxsds v3.0.4 scored;
libsais 2.9.1 measured, 2.12.1 scored).
