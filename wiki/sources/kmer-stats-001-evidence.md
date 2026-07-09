---
type: source
title: "Evidence: KMER-STATS-001 (K-mer statistics — composition summary + Shannon k-entropy)"
tags: [validation, analysis]
doc_path: docs/Evidence/KMER-STATS-001-Evidence.md
sources:
  - docs/Evidence/KMER-STATS-001-Evidence.md
source_commit: bb4c7f6095c6934658109faa87d4e9a6734ca72a
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: KMER-STATS-001

The validation-evidence artifact for test unit **KMER-STATS-001** — **k-mer statistics**
(`KmerAnalyzer.AnalyzeKmers`): the summary-statistics view over a sequence's k-mer count profile,
returning a `KmerStatistics` bundle of `TotalKmers`, `UniqueKmers` (= **distinct** count),
`MaxCount`, `MinCount`, `AverageCount`, and **`Entropy`** (Shannon k-entropy in bits). This is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the operation itself is synthesized in [[k-mer-statistics]]. See [[test-unit-registry]] for how units
are tracked.

It is a **companion summary layer** over the shared k-mer counting primitive — `AnalyzeKmers`
consumes the same `CountKmers` multiset the rest of the K-mer family is built on and reduces it to a
few scalars, adding the one genuinely new formula in this unit: the **k-entropy** (Shannon entropy of
the k-mer frequency distribution).

## What this file records

- **Online sources:**
  - **Wikipedia — K-mer** (rank 4) — the **total** count `L − k + 1` (overlapping length-k windows)
    and universe size `n^k`; the `AGAT` and `GTAGAGCTGT` worked count tables (k=1/2/3 totals
    10/9/8; distinct 4/7/8) that pin `TotalKmers`, `UniqueKmers`, `MaxCount`, `MinCount`.
  - **BioInfoLogics — k-mer counting part I (Clavijo 2018)** (rank 4) — corroborates `L − k + 1` and
    the **distinct vs unique** distinction; the `ATCGATCAC` k=3 table (ATC=2, others=1 → distinct 6,
    singletons 5). Used here only for the distinct count (the singleton/"unique" count is a *separate*
    unit, KMER-UNIQUE-001).
  - **Spectral concepts in genome informational analysis (Manca et al. 2021, arXiv:2106.15351)**
    (rank 1) — the **k-entropy** definition: `E_k(G) = −Σ p(α) log₂ p(α)` with
    `p(α) = mult(α) / (|G| − k + 1)` — the relative frequency of each k-mer over the `L − k + 1`
    windows.
  - **Entropy–Rank Ratio (arXiv:2511.05300)** (rank 1) — corroborates the single-sequence form
    `H_k(s) = −Σ pᵢ log₂ pᵢ` (base-2 / bits, convention `0·log 0 = 0`).

- **Documented corner cases / failure modes:** `k > L` (window doesn't fit, `L − k + 1 ≤ 0`) → no
  k-mers → all-zero statistics; **homopolymer / single distinct k-mer** (e.g. `AAAA`, k=2 → `AA×3`)
  → single-component distribution `p=1` → entropy `−1·log₂1 = 0` (minimum diversity); **all k-mers
  distinct** → each `p = 1/(L−k+1)` → entropy `log₂(L−k+1)` (maximum diversity at fixed total).

- **Datasets (documented oracles):**
  - `GTAGAGCTGT` (L=10): **k=1** → total 10, distinct 4, max 4 (G), min 1 (C), average 10/4 = 2.5,
    entropy 1.846439344671… bits; **k=2** → 9 / 7 / 2 / 1, average 9/7 ≈ 1.29; **k=3** → 8 / 8 / 1 / 1,
    entropy `log₂8 = 3.0`.
  - `ATCGATCAC` (L=9), k=3 → total 7, distinct 6, max 2 (ATC), min 1, average 7/6 ≈ 1.17, entropy
    2.521640636343… bits.
  - `AAAA` (L=4), k=2 → total 3, distinct 1, max=min=3, average 3.0, entropy 0 (homopolymer corner).

- **Test-coverage recommendations:** MUST — `TotalKmers = L − k + 1`; `UniqueKmers` = distinct count;
  `MaxCount`/`MinCount` = observed multiplicity extremes; `AverageCount = total/distinct`;
  `Entropy = −Σ p log₂ p`; invariants `TotalKmers == Σ counts` and `UniqueKmers == distinct(CountKmers)`
  cross-checked. SHOULD — homopolymer → entropy 0, max==min==total; all-distinct → entropy
  `log₂(distinct)`, max==min==1; empty / `k > L` → all-zero. COULD — case-insensitivity; `k ≤ 0` →
  `ArgumentOutOfRangeException`.

## Deviations and assumptions

Two **assumptions**, both presentation-only and non-correctness-affecting:

- **`AverageCount` rounded to 2 decimals** (`Math.Round(averageCount, 2)`): the literature defines
  average multiplicity as the exact ratio `total/distinct`; the display rounding is a repository
  contract choice (tests assert the rounded value; the exact ratio is verified per dataset row).
- **`Entropy` reported unrounded in bits** (log base 2, per the k-entropy sources): tests assert exact
  values within `1e-10`.

No source contradictions — Wikipedia/BioInfoLogics count tables and the Manca / Entropy–Rank Ratio
k-entropy formulas are mutually consistent.
