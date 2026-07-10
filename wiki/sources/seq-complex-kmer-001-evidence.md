---
type: source
title: "Evidence: SEQ-COMPLEX-KMER-001 (k-mer Shannon entropy — the complexity-family entry point)"
tags: [validation, analysis]
doc_path: docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md
sources:
  - docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md
source_commit: 49a4b6f93203c68a4ef386c8867dce3dcabf99d9
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-COMPLEX-KMER-001

The validation-evidence artifact for test unit **SEQ-COMPLEX-KMER-001** — **k-mer entropy**, the
Shannon entropy `H = −Σ pᵢ log₂ pᵢ` of the overlapping k-mer frequency distribution
(`pᵢ = nᵢ/(L−k+1)`), in **bits**. It is the **entropy member of the `SEQ-COMPLEX-*` sequence
complexity family** (siblings [[sequence-complexity-compression-lempel-ziv]] LZ76 and
[[dust-low-complexity-score]] DUST), implemented as `SequenceComplexity.CalculateKmerEntropy`.
The measure is **mathematically identical** to the Shannon k-entropy already written up as the
`Entropy` field of [[k-mer-statistics]] — so this file records the source trace and worked oracles,
and the method is enriched onto that concept page rather than duplicated. One instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; see
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Li, H. (2025)** *Finding low-complexity DNA sequences with longdust* arXiv:2509.07357 (rank 1,
    the same primary behind the DUST unit) — the operational definition used here: k-mers are the
    `L − k + 1` **overlapping** sliding-window words; `H = −Σ pᵢ log₂ pᵢ` with `pᵢ = nᵢ/N`,
    `N = L − k + 1`; low-complexity ⇒ skewed distribution ⇒ **low** entropy, high-complexity ⇒
    uniform ⇒ **high** entropy.
  - **Çakır et al. (2025)** *Entropy–Rank Ratio* arXiv:2511.05300 (rank 1, also corroborates
    [[k-mer-statistics]]) — the same entropy formula `S(w) = −Σ bⱼ log bⱼ`, `bⱼ = aⱼ/M`,
    **base-2 → bits**, single-nucleotide max = `log₂4 = 2` bits; saturation to `log(λ)` for long
    uniform i.i.d. input.
  - **Shannon (1948)** via secondary expositions (rank 4) — the well-established **bounds**
    `0 ≤ H ≤ log(k)`: `H = 0` for a deterministic (certain) distribution, `H = log_b(n)` for a
    uniform distribution over `n` outcomes.
- **Datasets (hand-derived worked oracles, `N = L−k+1`, bits):**

  | Input | k | Counts | N | H (bits) |
  |-------|---|--------|---|----------|
  | `ACGT` | 1 | A,C,G,T each 1 | 4 | `log₂4 = 2.0` (uniform) |
  | `ACGT` | 2 | AC,CG,GT each 1 | 3 | `log₂3 = 1.5849625007211562` (all-distinct) |
  | `ATATAT` | 2 | AT=3, TA=2 | 5 | `0.9709505944546686` (binary entropy of p=0.6) |
  | `AAAA` | 2 | AA=3 | 3 | `0.0` (deterministic / homopolymer) |
  | `AAACGT` | 2 | AA=2, AC=1, CG=1, GT=1 | 5 | `1.9219280948873623` (= log₂5 − 0.4) |
  | `AC` | 5 | none (L < k) | — | `0.0` |

- **Documented corner cases:** single repeated k-mer (homopolymer / tandem repeat) ⇒ `p = 1` ⇒
  `H = 0` (lowest complexity); all-distinct k-mers ⇒ `pᵢ = 1/N` ⇒ maximum `H = log₂N`.
- **Recommended coverage:** the five worked entropy oracles exactly; `L < k → 0`; invalid k / null
  guards; string-vs-`DnaSequence` overload agreement (case-insensitive); invariant
  `0 ≤ H ≤ log₂(L−k+1)` (Shannon bounds).

## Deviations and assumptions

Two **ASSUMPTIONs**, both API-contract shape only — neither changes an entropy value for valid input:

1. **`L < k` (sequence shorter than k) returns 0** — no source numerically specifies it (no k-mers
   exist; entropy of an empty multiset is conventionally 0). Resolved to match the `SequenceComplexity`
   siblings (`CalculateLinguisticComplexity`, `EstimateCompressionRatio`), which return 0 for
   empty/too-short input.
2. **Invalid `k (< 1)` → `ArgumentOutOfRangeException`; null `DnaSequence` → `ArgumentNullException`;
   null/empty string → 0** — library-API contract, not entropy literature; matched to the sibling
   method guards in the same class.

The implementation (`SequenceComplexity.CalculateKmerEntropyCore`) is a verbatim realization of
`H = −Σ (nᵢ/N) log₂(nᵢ/N)` over the overlapping k-mer multiset — the same formula as
`KmerAnalyzer.AnalyzeKmers.Entropy`.

## Contradictions

**None.** Li 2025, the Entropy–Rank Ratio preprint, and the Shannon expositions agree on the formula,
the base-2 (bits) convention, and the `0 ≤ H ≤ log₂N` bounds. This unit and [[k-mer-statistics]]
compute the identical measure from two different classes — not a contradiction, a duplicate entry
point (noted on that concept page).
