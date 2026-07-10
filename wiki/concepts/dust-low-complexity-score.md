---
type: concept
title: "DUST low-complexity score (triplet-frequency masking score)"
tags: [analysis, algorithm]
mcp_tools:
  - complexity_dust_score
  - complexity_mask_low
  - dust_score
  - mask_low_complexity
sources:
  - docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md
source_commit: dfe172057b03d328680fbf256469d7675c2604a4
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-complex-dust-001-evidence
      evidence: "Test Unit ID: SEQ-COMPLEX-DUST-001 ... Algorithm: DUST Score (triplet-frequency low-complexity score of Morgulis et al. 2006 SDUST/DUST)"
      confidence: high
      status: current
---

# DUST low-complexity score (triplet-frequency masking score)

The **DUST** score of Morgulis et al. (2006) — the standard **low-complexity DNA masking** score,
computed from **nucleotide triplet frequencies**. For a string `x` of length `L`:

```
score(x) = ( ∑_t c_x(t)·(c_x(t)−1)/2 ) / (L − 2)
```

where `c_x(t)` is the number of times triplet `t` occurs in `x` and the sum runs over all triplets
(`k = 3` hardcoded; `L − 2` = the number of triplets in `x`). Repeated triplets raise
`∑ c(c−1)/2`, so a **HIGH score indicates LOW complexity** — the high-scoring regions are the ones
masked. Validated as test unit **SEQ-COMPLEX-DUST-001**; the source trace and worked oracles are in
[[seq-complex-dust-001-evidence]], [[test-unit-registry]] tracks the unit, and see
[[algorithm-validation-evidence]] for the artifact pattern.

## Where it sits in the complexity family

DUST is the **triplet-frequency** member of the sequence complexity/entropy family, and the
**DNA-sequence low-complexity masker** specifically. It is a *distinct scalar measure* from its
siblings, not a re-derivation:

- vs. **compression complexity** — [[sequence-complexity-compression-lempel-ziv]] (Lempel–Ziv LZ76
  phrase count) measures *adaptive variable-length pattern buildup along the whole sequence*; DUST
  instead sums over a **fixed `k = 3` triplet-count distribution** in a window. Both live in the same
  `SEQ-COMPLEX-*` family and both are low exactly where repeats are, but DUST is directly a
  **masking** score with a published threshold, LZ a general normalized complexity scalar.
- vs. **protein SEG** — [[protein-low-complexity-seg]] is the **protein counterpart** low-complexity
  detector, but over a 20-letter amino-acid alphabet and using **Shannon entropy of composition**
  (`H = −Σ pᵢ log₂ pᵢ`, *low* H ⇒ low complexity). DUST is the DNA analogue with the opposite
  numeric direction (*high* DUST score ⇒ low complexity) and a triplet-frequency rather than
  single-residue-frequency statistic.
- vs. **k-mer k-entropy** — the family's **entropy member** ([[k-mer-statistics]], validated
  standalone as SEQ-COMPLEX-KMER-001 = `SequenceComplexity.CalculateKmerEntropy`) reduces a k-mer
  count profile to a Shannon k-entropy; DUST uses the same idea of a k-mer (triplet) count profile
  but a **sum-of-pair-counts** `∑ c(c−1)/2` statistic tuned for masking rather than an entropy, and
  in the **opposite numeric direction** (high DUST ⇒ low complexity, whereas low entropy ⇒ low
  complexity).
- vs. **explicit repeats** — [[repetitive-element-detection]] finds and types explicit repeated
  substrings; DUST collapses low-complexity/repetitiveness into a single maskable scalar. Link
  DNA low-complexity masking there as the repeats-family anchor; DUST is the score behind it.
- vs. **windowed complexity profile** — [[windowed-sequence-complexity-profile]] is the family's
  *profiling* member (SEQ-COMPLEX-WINDOW-001): a sliding-window scan emitting a per-position
  `ComplexityPoint` profile (Shannon entropy + linguistic complexity per window) that a caller
  thresholds to locate low-complexity regions — the *scanning* companion to DUST's masking scalar,
  but built on Shannon + linguistic-complexity statistics rather than the `∑ c(c−1)/2` triplet sum.

## Default parameters and the masking threshold

| Parameter | Value | Source |
|-----------|-------|--------|
| Triplet / word size `k` | **3** | Morgulis 2006; Li 2025 (`k = 3` hardcoded) |
| Default window size `w` | **64** | Morgulis 2006 ("64-base windows"); lh3/sdust `int W = 64;` |
| Default threshold (mask if score >) | **2.0** | Li 2025 (level 20, score 2.0); lh3/sdust `T = 20` |

The reference implementation compares `if (rw * 10 > L * T)` — i.e. `rw/L > T/10` — so with `T = 20`
the mask test is exactly `score > 2.0`. The 2006 rewrite (SDUST) kept this scoring function
unchanged and only made the **masking rule symmetric and context-insensitive**. The internal
scoring form `S_S = (1/ℓ)·∑_t c(c−1)/2 − T` subtracts the threshold `T` from the raw score.

The reference incremental accumulation (`rw += cw[t]++` before increment) sums the pre-increment
counts `0+1+…+(c−1) = c(c−1)/2` per triplet, so it is provably equal to the closed-form sum above.

## Worked oracles (k = 3, divisor = L − 2)

| Input | L | Σ c(c−1)/2 | L−2 | Score |
|-------|---|-----------|-----|-------|
| `ATGC` | 4 | 0 (ATG=1, TGC=1) | 2 | **0.0** (all-distinct triplets ⇒ max complexity) |
| `ACGTACGT` | 8 | 2 (ACG=2, CGT=2) | 6 | **0.333…** |
| `AAAAAA` | 6 | 6 (AAA=4) | 4 | **1.5** |
| `ACACACAC` | 8 | 6 (ACA=3, CAC=3) | 6 | **1.0** |
| `AAAAAAAAAA` | 10 | 28 (AAA=8) | 8 | **3.5** |

## Corner cases and repository generalizations

- **Below triplet length (`L < 3`):** the score is **undefined** (`L−2 ≤ 0`, no triplets exist).
- **All-distinct triplets:** every `c(c−1)/2 = 0` ⇒ score **0** (maximum complexity).
- **Maximally repetitive homopolymer** (length `L`): one triplet repeated `L−2` times ⇒ maximal
  score `(L−2)(L−3)/2 / (L−2) = (L−3)/2`.
- **General word size (ASSUMPTION):** the repository exposes a `wordSize` parameter; for
  `wordSize = w` the normalization generalizes to the number of words `L − w + 1` (= `L − 2` at
  `w = 3`). Only `k = 3` is source-backed; exact oracle values are asserted for k = 3 only.
- **Input shorter than one word (`L < wordSize`, ASSUMPTION):** returns **0** (no repeats ⇒ minimal
  complexity) — a defined-output convention, not a source value.
- API-contract behaviours: case-insensitive (`ToUpperInvariant`); null DnaSequence →
  ArgumentNullException; null/empty string → 0; DnaSequence and string overloads agree.

## References

Morgulis A., Gertz E.M., Schäffer A.A. & Agarwala R. (2006) *J Comput Biol* 13(5):1028–1040
(DOI 10.1089/cmb.2006.13.1028, PMID 16796549); Li H. (2025) *Finding low-complexity DNA sequences
with longdust* arXiv:2509.07357; reference implementation lh3/sdust `sdust.c`. Full citations in
[[seq-complex-dust-001-evidence]] (do not duplicate here). A
[[research-grade-limitations|research-grade]] implementation of the standard DUST score.
</content>
