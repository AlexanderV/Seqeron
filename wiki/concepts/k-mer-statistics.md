---
type: concept
title: "K-mer statistics (composition summary + Shannon k-entropy over the count profile)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/KMER-STATS-001-Evidence.md
  - docs/algorithms/K-mer/K-mer_Statistics.md
source_commit: bb4c7f6095c6934658109faa87d4e9a6734ca72a
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-stats-001-evidence
      evidence: "Test Unit ID: KMER-STATS-001; Algorithm: K-mer Statistics (docs/algorithms/K-mer/K-mer_Statistics.md)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:asynchronous-kmer-counting
      source: kmer-stats-001-evidence
      evidence: "AnalyzeKmers is a summary layer over the same synchronous CountKmers sliding-window count multiset the K-mer family shares (KmerAnalyzer); UniqueKmers == distinct count from CountKmers, TotalKmers == sum of counts == L-k+1"
      confidence: high
      status: current
---

# K-mer statistics (composition summary + Shannon k-entropy over the count profile)

The **K-mer** family's *summary-statistics* operation: `KmerAnalyzer.AnalyzeKmers` reduces a
sequence's k-mer count profile to a small `KmerStatistics` bundle — `TotalKmers`, `UniqueKmers`,
`MaxCount`, `MinCount`, `AverageCount`, and `Entropy` (Shannon **k-entropy** in bits). It is a
**companion summary layer**, not a new counting algorithm: it calls the same `CountKmers`
sliding-window primitive the rest of the family shares and aggregates the resulting multiset. Its one
genuinely distinct piece of content is the **k-entropy** — the Shannon entropy of the k-mer frequency
distribution, a diversity/complexity measure with its own literature. Validated under test unit
**KMER-STATS-001** (record [[kmer-stats-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Core model

For a sequence `S` of length `L` and window length `k` (`1 ≤ k ≤ L`), the k-mers are the multiset of
the `L − k + 1` overlapping length-`k` windows. Let `mult(α)` be the multiplicity of k-mer `α` and
`D` the number of **distinct** k-mers. `AnalyzeKmers` computes:

| Field | Definition | Notes |
|-------|------------|-------|
| `TotalKmers` | `L − k + 1` = `Σ mult(α)` | number of overlapping windows (Wikipedia) |
| `UniqueKmers` | `D` = distinct k-mer count | **the distinct count**, *not* singletons — see the gotcha below |
| `MaxCount` | `max_α mult(α)` | most-frequent k-mer's multiplicity |
| `MinCount` | `min_α mult(α)` | least-frequent k-mer's multiplicity |
| `AverageCount` | `TotalKmers / D` = `(L − k + 1)/D` | average multiplicity; **rounded to 2 dp** in the output |
| `Entropy` | `−Σ p(α) log₂ p(α)`, `p(α) = mult(α)/(L − k + 1)` | Shannon **k-entropy**, bits, unrounded |

The **k-entropy** is the one formula unique to this unit: `E_k = −Σ p(α) log₂ p(α)` over the relative
frequencies of the k-mers (Manca et al. 2021, *Spectral concepts in genome informational analysis*;
corroborated by the Entropy–Rank Ratio preprint, arXiv:2511.05300), base-2 (bits), convention
`0·log 0 = 0`. It measures the diversity/evenness of the k-mer composition: minimal (0) for a single
distinct k-mer, maximal (`log₂(L−k+1)`) when every window is distinct.

## Gotcha: `UniqueKmers` is the *distinct* count, not the singleton count

Despite the name, the `UniqueKmers` field holds `D` = the number of **distinct** k-mers (each different
k-mer counted once), **not** the number of k-mers that occur exactly once. The "unique" (count == 1)
singleton notion is a **separate** unit, KMER-UNIQUE-001 — synthesized in
[[unique-and-mincount-kmers]] (`KmerAnalyzer.FindUniqueKmers`). For `ATCGATCAC` k=3 the distinct count
is 6 (`UniqueKmers = 6`) while the singleton count is 5.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `TotalKmers == Σ_α mult(α) == L − k + 1` (for `1 ≤ k ≤ L`) | one k-mer per start position |
| INV-02 | `UniqueKmers == D`, the distinct count from `CountKmers` | direct definition; cross-checked independently |
| INV-03 | `MinCount ≤ AverageCount·(rounding) ≤ MaxCount`, and `MaxCount == MinCount` ⇔ every distinct k-mer equally frequent | extremes bound the mean of the multiplicity distribution |
| INV-04 | single distinct k-mer (homopolymer) ⇒ `Entropy = 0`, `MaxCount == MinCount == TotalKmers` | one-component distribution `p = 1`, `−1·log₂1 = 0` |
| INV-05 | all k-mers distinct ⇒ `Entropy = log₂ D`, `MaxCount == MinCount == 1` | uniform distribution `p = 1/D` over `D` equiprobable k-mers |
| INV-06 | `k > L` or empty `S` ⇒ all-zero `KmerStatistics` | `L − k + 1 ≤ 0` ⇒ no k-mers |

## Contract and edge cases

- **Inputs:** `sequence` (string; empty/null → all-zero statistics; upper-cased internally so the
  statistics are **case-insensitive**), `k` (`int`, `> 0`; `k > L` → all-zero result). `k ≤ 0` →
  `ArgumentOutOfRangeException` (shared with the sibling counting methods).
- **Output:** a `KmerStatistics` value with the six fields above; `AverageCount` rounded to 2 decimals,
  `Entropy` unrounded in bits.
- **Edge cases:** empty / `k > L` → all fields 0; homopolymer → entropy 0 (INV-04); all-distinct →
  entropy `log₂ D` (INV-05); lower/mixed case → same statistics as upper.

## Worked oracles

- `GTAGAGCTGT` (L=10) **k=1** → Total 10, Unique 4, Max 4 (G), Min 1 (C), Average 2.5, Entropy
  1.846439344671… bits.
- `GTAGAGCTGT` **k=2** → 9 / 7 / 2 / 1, Average 9/7 ≈ 1.29 (rounded). **k=3** → 8 / 8 / 1 / 1, Entropy
  `log₂8 = 3.0` (8 equiprobable k-mers, INV-05).
- `ATCGATCAC` (L=9) k=3 → Total 7, Unique 6, Max 2 (ATC), Min 1, Average 7/6 ≈ 1.17, Entropy
  2.521640636343… bits.
- `AAAA` (L=4) k=2 → Total 3, Unique 1, Max=Min=3, Average 3.0, Entropy 0 (homopolymer, INV-04).

## Deviations and assumptions

- **No source contradictions; deviations = None.** The Wikipedia/BioInfoLogics count tables and the
  Manca / Entropy–Rank Ratio k-entropy formula are implemented verbatim.
- **ASSUMPTION (presentation-only):** `AverageCount` is rounded to 2 decimals (`Math.Round(_, 2)`) —
  the underlying `total/distinct` ratio is exact; `Entropy` is returned unrounded in bits (log base 2).
  Neither affects correctness; tests assert the rounded average and entropy within `1e-10`.

## Relation to other units

- **Summary layer over the shared count.** `AnalyzeKmers` aggregates the same synchronous
  `KmerAnalyzer.CountKmers` sliding-window multiset that underlies [[asynchronous-kmer-counting]] and
  [[both-strand-kmer-counting]] (a future sync-count concept will own the `L − k + 1` count
  definition, and this page will link to it) — it computes *statistics over* the profile rather than
  producing it.
- Distinct from the family's other operations: [[k-mer-generation]] enumerates the full `n^k`
  **universe** (sequence-independent), [[k-mer-positions]] indexes **where** one k-mer occurs, and
  [[k-mer-euclidean-distance]] / [[kmer-jaccard-similarity]] compare the frequency/presence profiles of
  **two** sequences — whereas this unit summarizes the profile of **one**.
- The **k-entropy** here is the k-mer specialization of Shannon sequence entropy; a general
  `shannon-entropy` unit (pending ingest) would relate to it as the character-level counterpart.
  A **different** member of the same complexity/entropy family is the compression-based
  [[sequence-complexity-compression-lempel-ziv]] (Lempel–Ziv LZ76 count): it counts adaptively
  discovered **variable-length** phrases along the whole sequence rather than a fixed-`k`
  frequency distribution, so it captures ordered pattern buildup that a fixed-`k` k-entropy misses.
