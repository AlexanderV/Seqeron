---
type: source
title: "Evidence: ASSEMBLY-STATS-001 (Assembly statistics — N50/L50/auN)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-STATS-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-STATS-001-Evidence.md
source_commit: ca833349c871b6ae6b0fc30b5a251408d2ede693
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-STATS-001

The validation-evidence artifact for test unit **ASSEMBLY-STATS-001** — assembly statistics:
the contiguity/QC summary metrics (N50, L50, Nx, Lx, N90, L90, auN) plus totals, largest/smallest,
GC, and an N-run gap summary computed over a set of contig lengths. One instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the metric definitions,
the inclusive-≥ threshold, the auN formula, the published oracles and the two assumptions are
summarized in [[assembly-statistics]], the anchor for the assembly STATS family. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (all accessed 2026-06-13):
  - **Miller, Koren & Sutton (2010) — "Assembly algorithms for next-generation sequencing data"**
    (rank 1; *Genomics* 95(6):315–327, via WebFetch of the PMC full text) — the canonical N50
    definition (§1.2): the contig N50 is "the length of the smallest contig in the set that contains
    the fewest (largest) contigs whose combined length represents at least 50% of the assembly" →
    largest-first, **cumulative ≥ 50%** (inclusive) boundary.
  - **Wikipedia — "N50, L50, and related statistics"** (rank 4; rendered article + `action=raw`
    wikitext) — used for its worked examples and as a cross-check of the Miller primary: N50 = a
    **length** ("shortest contig at 50% of total assembly length" / weighted median), L50 = a
    **count**, N90 defined analogously at 90% (N90 ≤ N50), and the Assembly A / Assembly B worked
    examples.
  - **QUAST reference implementation — `quast_libs/N50.py`** (rank 3; ablab/quast master, raw source
    read) — the cumulative comparison `NG50_and_LG50` (`limit = ref·(100−pct)/100`, stop when
    `s <= limit`, i.e. cumulative ≥ threshold; `N50()` = `NG50(numlist, sum(numlist))`, `L50()` =
    iteration count) and `au_metric` = `sum(n² for n in numlist) / denum` → **auN = ΣLᵢ²/ΣLⱼ**
    (empty list → `None`).
  - **Heng Li (2020) — "auN: a new metric to measure assembly contiguity"** (rank 3; the originating
    definition, samtools/minimap2 author) — Nx: "Contigs no shorter than Nx cover x% of the
    assembly"; **auN = ∑ᵢ Lᵢ·(Lᵢ/∑ⱼ Lⱼ) = ∑ᵢ Lᵢ²/∑ⱼ Lⱼ**, the area under the Nx curve, more stable
    than the discrete N50.
- **Datasets (published oracles):**
  - **Assembly A** {80,70,50,40,30,20}, total 290 → **N50 = 70, L50 = 2**, N90 = 30, L90 = 5,
    auN = 16700/290 ≈ 57.586 (Wikipedia worked example, cross-checked vs QUAST).
  - **Assembly B** {80,70,50,40,30,20,10,5}, total 305 → **N50 = 50, L50 = 3** (threshold-shift
    example, same Wikipedia article).
  - **auN formula check** {100,80,60,40,20}, total 300 → auN = 22000/300 = 73.333…, N50 = 80, L50 = 2
    (QUAST `au_metric` + lh3 2020).
- **Corner cases / failure modes** — inclusive threshold boundary (N50 is the *first* largest-first
  contig at which cumulative length reaches **at least** 50%, confirmed by Miller's "at least 50%"
  and QUAST's `s <= limit`); monotonicity **N90 ≤ N50**, **L90 ≥ L50**; empty input has no defined
  N50/L50/auN (QUAST returns `None`, the repository returns 0).
- **Recommended coverage** — MUST: N50=70/L50=2 on Assembly A; N50=50/L50=3 on Assembly B; the
  inclusive boundary (cumulative reaching exactly 50% triggers Nx); auN = ΣL²/ΣL exact (73.333… on
  {100,80,60,40,20}); N90 ≤ N50 and L90 ≥ L50; `FindGaps` detects N-runs with exact 0-based
  inclusive `[Start,End]` + length. SHOULD: `CalculateStatistics` aggregates the full summary on a
  multi-contig assembly; `FindGaps` `minGapLength` filter + leading/trailing gaps. COULD: empty input
  returns zeros without throwing.

## Assumptions (from the artifact)

Two assumption records, both outside the source-backed N50/L50/Nx/auN contract:

1. **Empty-input return value (non-correctness-affecting).** QUAST returns `None` for an empty
   contig set (`au_metric` asserts `len > 0`; `NG50_and_LG50` returns `(None, None)`); the repository
   returns the all-zero `AssemblyStatistics` with `Nx = Lx = 0` / `auN = 0`. An API-shape choice (no
   valid N50 exists for an empty assembly), documented in the algorithm doc §6.1 — it changes no
   defined value.
2. **Median convention in `CalculateStatistics.MedianLength`.** The N50 sources define no assembly
   "median contig length"; the repository reports the **upper median** (`lengths[count/2]` over the
   descending-sorted list). This field is auxiliary and outside the cited contract — tested as
   implemented behavior, flagged as not source-derived. The canonical statistics (N50, L50, N90, L90,
   largest, smallest, totals, GC, gaps) are source-backed.

No contradictions among the sources — Miller 2010, the Wikipedia worked examples, QUAST `N50.py`, and
Heng Li (2020) give the identical largest-first, inclusive-≥ cumulative-threshold definitions, and
QUAST's `au_metric` matches Heng Li's ΣL²/ΣL exactly.
