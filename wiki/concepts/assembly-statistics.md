---
type: concept
title: "Assembly statistics (N50 / L50 / Nx / Lx / auN)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-STATS-001-Evidence.md
  - docs/algorithms/Assembly/Assembly_Statistics.md
  - docs/Validation/reports/ASSEMBLY-STATS-001.md
source_commit: d584af4da843a888434b5c54e7277e8f3085b085
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-stats-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-STATS-001 ... Assembly Statistics (N50 / L50 / Nx / Lx / auN, gap detection, contiguity summary)"
      confidence: high
      status: current
---

# Assembly statistics (N50 / L50 / Nx / Lx / auN)

**Assembly statistics** are the contiguity/QC summary metrics computed over a finished set of
contig (or scaffold) lengths — the report you read *after* assembly to judge how contiguous the
result is. Unlike the contig-*building* steps ([[de-bruijn-graph-assembly|DBG]] /
[[overlap-layout-consensus-assembly|OLC]] / [[contig-merge-overlap-collapse|merge]]) and the
downstream [[scaffolding]] layout step, this unit consumes lengths and emits numbers: **N50, L50,
Nx, Lx, N90, L90, auN**, plus totals, largest/smallest, GC, and an N-run **gap** summary. It is the
anchor for the assembly **STATS** family, validated under test unit **ASSEMBLY-STATS-001**. The
literature-traced validation record is [[assembly-stats-001-evidence]]; the independent two-stage
re-validation verdict (Stage A/B PASS, CLEAN) is [[assembly-stats-001-report]]; [[test-unit-registry]]
tracks the unit and [[algorithm-validation-evidence]] describes the artifact pattern.

## The metric definitions (source-traced)

Sort contig lengths **largest-first**; let `T = Σ Lᵢ` be the total assembly length.

- **Nx** — the length of the shortest contig such that contigs of that length or longer cover at
  least **x%** of `T`. Miller, Koren & Sutton (2010): the contig N50 is "the length of the smallest
  contig in the set that contains the fewest (largest) contigs whose combined length represents at
  least 50% of the assembly." Heng Li (2020): "Contigs no shorter than Nx cover x% of the assembly."
- **Lx** — the **count** of contigs needed to reach that x% threshold (Wikipedia: "count of smallest
  number of contigs whose length sum makes up half of genome size"). **Nx is a length; Lx is a
  count.** N50/L50 are the x=50 case; N90/L90 the x=90 case.
- **Inclusive (≥) threshold.** The cumulative sum reaches the threshold with **"at least" x%** — the
  boundary is inclusive. Confirmed independently by Miller's "at least 50%" and the QUAST reference
  stop test `s <= limit` (where `s = total − cumulative`, `limit = total·(100−x)/100`), i.e.
  cumulative ≥ threshold.
- **Monotonicity.** Raising the threshold cannot raise Nx: **N90 ≤ N50** and **L90 ≥ L50**.
- **auN** — the **area under the Nx curve**, a single continuous contiguity number that is more
  stable than the discrete N50. Heng Li (2020) / QUAST `au_metric`:
  **auN = Σᵢ Lᵢ·(Lᵢ / Σⱼ Lⱼ) = Σᵢ Lᵢ² / Σⱼ Lⱼ** (each length weighted by its own share of the total).

The auxiliary fields (largest, smallest, total length, contig count, GC%, and the gap summary) are
plain aggregates. **Gaps** are maximal runs of `N`/`n`, each reported as a 0-based inclusive
`[Start, End]` span plus length, subject to a `minGapLength` filter.

## Worked oracles (published)

Traced to the Wikipedia "N50, L50, and related statistics" worked examples (which cite Miller 2010)
and cross-checked against QUAST `N50()`/`L50()`/`au_metric` and Heng Li's auN formula:

| Dataset | Lengths | Total | N50 | L50 | N90 | L90 | auN |
|---------|---------|-------|-----|-----|-----|-----|-----|
| Assembly A | {80, 70, 50, 40, 30, 20} | 290 | **70** (80+70 > 145) | **2** | 30 (cum 270 ≥ 261) | 5 | 16700/290 ≈ **57.586** |
| Assembly B | {80, 70, 50, 40, 30, 20, 10, 5} | 305 | **50** (80+70+50 > 152.5) | **3** | — | — | — |
| auN check | {100, 80, 60, 40, 20} | 300 | 80 (100+80 ≥ 150) | 2 | — | — | 22000/300 = **73.333…** |

## Edge cases and assumptions (from the artifact)

- **Empty input — API-shape choice, non-correctness-affecting.** No N50/L50/auN is *defined* for an
  empty contig set: QUAST returns `None` (`au_metric` asserts `len > 0`; `NG50_and_LG50` returns
  `(None, None)`). The repository instead returns the all-zero `AssemblyStatistics` with
  `Nx = Lx = 0` / `auN = 0` rather than throwing. This changes no defined value; it is documented in
  the algorithm doc §6.1.
- **"Median contig length" is not source-derived.** The cited N50 literature defines no assembly
  "median contig length"; `CalculateStatistics.MedianLength` reports the **upper median**
  (`lengths[count/2]` over the descending-sorted list). This auxiliary field sits outside the cited
  N50/L50/Nx/auN contract and is tested as implemented behavior, flagged as not source-backed. The
  canonical statistics (N50, L50, N90, L90, largest, smallest, totals, GC, gaps) *are* source-backed.

No contradictions among the sources — Miller 2010, the Wikipedia worked examples, QUAST `N50.py`, and
Heng Li (2020) give the identical largest-first, inclusive-≥, cumulative-threshold definitions;
QUAST's `s <= limit` and Miller's "at least 50%" are the same inclusive boundary, and QUAST's
`au_metric` matches Heng Li's ΣL²/ΣL exactly.

## Relation to the other assembly steps

Assembly statistics is a **read-only summary** over finished lengths: it neither builds nor lays out
sequence. It sits **downstream** of [[de-bruijn-graph-assembly|DBG]] /
[[overlap-layout-consensus-assembly|OLC]] / [[contig-merge-overlap-collapse|merge]] contig
construction and is orthogonal to [[scaffolding]] (which introduces the `N`-gaps this unit then
*reports* on) and to per-base [[coverage-depth-calculation|coverage]] (a depth profile, not a
length-distribution metric). N50/auN are the headline numbers a coverage- or scaffold-improving
step is trying to move.
