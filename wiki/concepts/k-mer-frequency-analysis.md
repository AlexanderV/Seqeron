---
type: concept
title: "K-mer frequency analysis (normalized frequencies + count-of-counts spectrum + k-entropy over the shared count)"
tags: [analysis, algorithm]
sources:
  - docs/algorithms/K-mer/K-mer_Frequency_Analysis.md
source_commit: 6b60958bf165801aa5dd9c92ffed92970bac80ab
created: 2026-07-13
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-freq-001-report
      evidence: "Test Unit ID: KMER-FREQ-001; Algorithm: K-mer Frequency Analysis (docs/algorithms/K-mer/K-mer_Frequency_Analysis.md), KmerAnalyzer.GetKmerFrequencies / GetKmerSpectrum / CalculateKmerEntropy"
      confidence: high
      status: current
---

# K-mer frequency analysis (normalized frequencies + count-of-counts spectrum + k-entropy over the shared count)

The **distribution-shaping** layer of the K-mer family: three `KmerAnalyzer` methods that turn the raw
count multiset produced by [[k-mer-counting]] (`CountKmers`) into the three standard views of a k-mer
*distribution* — the **normalized frequency vector**, the **k-mer spectrum** (count-of-counts
histogram), and the **Shannon k-entropy** of the distribution. None of the three re-counts anything:
all delegate to `CountKmers`, so the shared `L − k + 1` count definition and its contract are owned by
[[k-mer-counting]]. Validated under test unit
[[kmer-freq-001-report|KMER-FREQ-001]]; [[test-unit-registry]] tracks the unit
and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Core model

For a sequence `S` of length `L`, window length `k` (`1 ≤ k ≤ L`), and the observed k-mer counts
`c_i = Count(w_i)` over the `L − k + 1` overlapping windows:

| Output | Method | Definition |
|--------|--------|------------|
| **Normalized frequencies** | `GetKmerFrequencies(string, int)` → `Dictionary<string,double>` | `f_i = c_i / Σ_j c_j` — each observed k-mer's relative frequency in `[0, 1]` (divided by the total observed count `Σ_j c_j = L − k + 1`, **not** by the `nᵏ` theoretical k-mer space) |
| **K-mer spectrum** | `GetKmerSpectrum(string, int)` → `Dictionary<int,int>` | the **count-of-counts** histogram: `count → number of distinct k-mers observed with that count` (the count dictionary inverted) |
| **K-mer entropy** | `CalculateKmerEntropy(string, int)` → `double` | Shannon entropy in **bits**: `H = −Σ_i f_i log₂ f_i`, convention `0·log 0 = 0` (implemented by iterating only over observed non-zero frequencies) |

Sources: Wikipedia (K-mer; Entropy), Shannon 1948, Rosalind (K-mer Composition), Teeling et al. 2004
(TETRA tetranucleotide signatures), Chor et al. 2009 (genomic k-mer spectra).

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Σ_i f_i = 1.0` whenever at least one k-mer exists | each count is divided by the total observed count |
| INV-02 | spectrum total `Σ_c (c × mult(c)) = L − k + 1` for `1 ≤ k ≤ L` | the spectrum re-buckets the exact counts, whose sum is `L − k + 1` (INV-01 of [[k-mer-counting]]) |
| INV-03 | `0 ≤ H ≤ log₂(D)`, `D` = distinct k-mer count | `H` is the Shannon entropy of a discrete distribution over `D` outcomes: `0` for one distinct k-mer, `log₂ D` when all `D` are equiprobable |

## Contract and edge cases

- **Inputs:** `sequence` (`string`; null/empty → empty frequency map, empty spectrum, entropy `0.0`),
  `k` (`int`; `k > L` → empty/empty/`0.0`). All three methods delegate input handling to `CountKmers`,
  so `k ≤ 0` on a **non-empty** sequence throws `ArgumentOutOfRangeException` through the underlying
  count routine (the same `k ≤ 0` validation-order behavior documented on [[k-mer-counting]]).
  Counting is case-insensitive (uppercased internally by `CountKmers`).
- **Outputs:** `frequencies` (`Dictionary<string,double>`), `spectrum` (`Dictionary<int,int>`),
  `entropy` (`double`, bits — `Math.Log2`, **unrounded**).
- **Edge cases:**

| Case | frequencies | spectrum | entropy | Rationale |
|------|-------------|----------|---------|-----------|
| Empty sequence | `{}` | `{}` | `0.0` | no k-mers exist |
| `k > L` | `{}` | `{}` | `0.0` | no valid windows |
| Single possible k-mer | `{ w: 1.0 }` | `{ 1: 1 }` | `0.0` | one-outcome distribution |
| Homopolymer `AAAA`, `k = 2` | `{ AA: 1.0 }` | `{ 3: 1 }` | `0.0` | all `L − k + 1 = 3` windows identical → one k-mer with count 3 |

## Complexity

All three are `O(n)` time (`n = L − k + 1`), `O(u)` space (`u` = distinct observed k-mers), on top of
the single `CountKmers` scan they share: `GetKmerFrequencies` divides each count by the total,
`GetKmerSpectrum` iterates the count *values* to build the histogram, and `CalculateKmerEntropy` sums
`−f log₂ f` over the frequencies. Memory grows with the number of unique observed k-mers, inherited
from the counting primitive (no distribution smoothing, no normalization against the `nᵏ` k-mer space).

## The k-entropy is the same Shannon k-entropy owned by [[k-mer-statistics]]

`CalculateKmerEntropy` here is a **third entry point** to the identical formula already synthesized on
[[k-mer-statistics]] — `H = −Σ p(α) log₂ p(α)`, `p(α) = mult(α)/(L − k + 1)`, bits — the same measure
carried by the `AnalyzeKmers.Entropy` field (KMER-STATS-001) and by
`SequenceComplexity.CalculateKmerEntropy` (SEQ-COMPLEX-KMER-001, the complexity-family member). It is
**not a different metric**, just a standalone `KmerAnalyzer` surface; the derivation, bounds
(INV-03 here = INV-04/INV-05 there), and worked oracles live on [[k-mer-statistics]] and are not
duplicated. `GetKmerFrequencies` and `GetKmerSpectrum` are the genuinely distinct content of this unit.

## Deviations and assumptions

- **No source contradictions.** Frequency normalization by total observed count, spectrum as a
  multiplicity histogram, and base-2 Shannon entropy over the observed distribution are implemented
  verbatim from the cited theory (Shannon 1948; Wikipedia K-mer/Entropy; Teeling 2004; Chor 2009).
- **Deviation (accepted):** the original document described **4-decimal rounding** of the entropy for
  numerical stability, but the current `CalculateKmerEntropy` returns the **raw `double` sum** with full
  floating-point precision (no explicit rounding step). Confirmed from the source; does not affect
  correctness.
- **Limitation:** frequencies are relative to *observed* k-mers only — the distribution is **not**
  smoothed and **not** normalized against the theoretical `nᵏ` k-mer space, so frequencies from
  sequences of different composition/length are comparable only as relative-frequency vectors, not as
  absolute probabilities over the full k-mer universe.

## Relation to other units

- **Built on the shared count.** All three methods reduce the same `KmerAnalyzer.CountKmers` multiset
  owned by [[k-mer-counting]] (`Σ counts = L − k + 1`); this page adds the *distribution* views on top
  of it rather than a new counting algorithm.
- **Frequency vector → distance.** The `GetKmerFrequencies` normalized vector (`count ÷ sum`) is the
  exact per-word frequency vector [[k-mer-euclidean-distance]] consumes to compute its L2 distance
  between two sequences; the `kmer_frequencies` MCP tool is surfaced through that concept.
- **Spectrum → error detection / assembly QC.** The `GetKmerSpectrum` count-of-counts histogram is the
  input shape for k-mer-spectrum error detection and assembly QC — the low-count (singleton) tail flags
  sequencing errors, the coverage peak estimates depth; see [[kmer-spectrum-error-correction]] (which
  owns the `kmer_spectrum` MCP tool).
- **Summary vs distribution.** [[k-mer-statistics]] reduces the same count profile to a small scalar
  `KmerStatistics` bundle (total/unique/max/min/average + entropy); this unit instead returns the full
  per-k-mer frequency map and the multiplicity histogram. The k-entropy is shared between them.
- **Applications** (from the spec): genome-assembly QC via k-mer spectra, metagenomics binning via
  tetranucleotide-frequency signatures (Teeling 2004), alignment-free sequence comparison via k-mer
  frequency profiles, and sequencing-error detection from low-frequency k-mers.
