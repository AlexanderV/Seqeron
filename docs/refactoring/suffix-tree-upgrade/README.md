# Suffix-tree upgrade: competitor analysis and plan (2026-11-03)

Owner's goal: make the Seqeron suffix tree better, on average, than every competitor. Method: take all competitors,
take all their capabilities, score each intersection in a two-dimensional table, then derive a plan whose projected
average beats all of them. **Status (2026-11-23): S0, S1a and S1b are done and measured** (an enhanced-suffix-array backend, an external-memory
builder, a persistent container, a concurrent reader; the independent review of S1, work package 1b.10, is still to come). **S2-S8 are not approved and wait for the owner.**
The plan stays a plan: `PLAN.md` is unchanged and its projected columns below are those of 2026-11-03. The Seqeron column of the matrix was re-scored on 2026-11-23
(`analysis/review/1b11_rescore.md`); the competitor columns were not touched.

| file | content |
|---|---|
| [COMPETITOR_MATRIX.md](COMPETITOR_MATRIX.md) | 12 competitors + Seqeron: 24 properties (0-4, weighted) and 81 functions (0-2); ranking; gap list (competitors of 2026-11-03, Seqeron re-scored 2026-11-23) |
| [PLAN.md](PLAN.md) | stages S0-S8: scope, oracle, numeric acceptance gates, projected cells, cost, risks, order, hand-off items, decisions |
| [BENCHMARK.md](BENCHMARK.md) | measurements on one machine: Seqeron (in-memory, persistent), Go, libsais, libdivsufsort, sdsl-lite, rust-bio, pydivsufsort, Python suffix-trees, BWA, Bowtie2 |
| `analysis/` | rubric (`SPEC.md`, `FUNCTIONS.md`), per-cell evidence TSVs, Seqeron audit, benchmark harness and raw rows, hand-off package, scripts |

## Headline

| | Seqeron now (re-scored 2026-11-23; 1.80 and 55.6% on 2026-11-03) | best competitor | after the recommended core (S0-S6, S8), projection of 2026-11-03 | after all stages, projection of 2026-11-03 |
|---|---|---|---|---|
| property mean (24, weighted, 0-4) | 2.28 (rank 4 of 13) | 2.48 (ropebwt3) | 2.80 | 3.02 |
| function coverage (81) | 72.8% (56.8% on the tree layer only; was 39.5%) | 52.5% (sdsl) | 93.8% | 95.1% |
| first place over 2000 random weightings | not recomputed (0% on 2026-11-03) | - | 99.5% | 100% |

Measured at 10 Mbp random DNA on 2026-11-23 (5 fresh processes per row, median [min-max], one thread unless stated; `ADR-004` section 5a "1b.9", rows `analysis/bench/results_1b9_*`;
x86-64 Linux, 4 vCPU): the enhanced suffix array (suffix array, LCP, child table) is built in memory in 1.138 s [1.118-1.174] (the public `CreateEsa` 0.841 s [0.832-0.907]; libsais in the same session 0.344 s,
pooled median of 15 runs; the Ukkonen tree of 2026-11-03: 14.6 s), takes 10.00 B/char in memory (Ukkonen: 169) and 9.25 B/char (core) or 13.25 B/char (with the Psi section) as a persistent container file;
the building process peaks at 15.4 B/char (153.8 MB). It counts 647,176 20-mers/s [611,792-824,474] (best tool of 2026-11-03: 747k) and locates 990,228 8-mers/s [876,565-1,627,933] (best: 1.14M; 2,173,239/s with the lazy `Locate`).
On four threads over one shared tree it counts 3,098,222 20-mers/s (3.96x its own one-thread figure of 783,122 on the same instrument) and locates 4,841,674 8-mers/s. At 100 Mbp the container is built by the external-memory pipeline in
17.88 s [15.88-19.43] on four threads with a peak anonymous memory (RssAnon) of 253.2 MB [247.1-256.1] (286.3 MB on the genome-like text; the gate is 320 MB) and 13.25 B/char on disk; every one of 10 chains passed its 11 checks.
**Not met:** the first `Open` of a container in a fresh process takes 15.24 ms [13.09-15.37] against the 10 ms target (8.5 ms with ahead-of-time compilation), and the file-backed tree on a 20-letter text counts at 0.934 of the in-memory tree on four threads.
**Not measured:** 3x10^9 symbols (the owner's run is pending), arm64, macOS, Windows, the cold-cache latency of queries. The container builder is not public API yet.
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
