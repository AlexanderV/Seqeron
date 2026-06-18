# DUST Score

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-DUST-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

The DUST score is a heuristic measure of nucleotide-sequence low complexity introduced for masking repetitive DNA in BLAST [1]. It counts how often each overlapping triplet (3-mer) recurs in a sequence and combines those counts into a single score: highly repetitive sequences (e.g. homopolymers, simple satellites) score high, while sequences with all-distinct triplets score 0. It is a deterministic heuristic, not a probabilistic model; this unit implements the per-sequence/per-window complexity score, the quantity DUST/SDUST thresholds to decide masking.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Low-complexity regions (homopolymers, tandem repeats, simple sequence repeats) produce spurious alignments and inflate database-search noise. DUST ("Decreasing Uncertainty in Sequence Tags") masks such regions by scoring fixed-size windows and masking windows whose score exceeds a threshold [1]. The score function is unchanged between the original and the symmetric (SDUST) implementation; only the masking rule differs [1].

### 2.2 Core Model

For a sequence `x` of length `L`, let `c_t` be the number of occurrences of triplet `t` among the `L − 2` overlapping triplets. The complexity score is [2]:

```
S(x) = ( Σ_t  c_t·(c_t − 1)/2 ) / (L − 2)
```

The numerator `Σ_t c_t(c_t−1)/2` counts pairs of identical triplets; the denominator `L − 2` is the number of triplets [2]. Equivalently, the reference implementation accumulates the score incrementally as `rw += cw[t]++` when adding each triplet, so the running sum equals `Σ_t c_t(c_t−1)/2` [3]. The thresholded form is `S_S(x) = (1/ℓ) Σ_t c_t(c_t−1)/2 − T` with ℓ = L − 2 [2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `S(x) ≥ 0` | Each `c(c−1)/2 ≥ 0` and the divisor `L − 2 > 0` for `L ≥ 3` [2]. |
| INV-02 | All-distinct triplets ⇒ `S(x) = 0` | Every `c_t = 1` ⇒ `c(c−1)/2 = 0` [2]. |
| INV-03 | Homopolymer of length `L` ⇒ `S = (L−3)/2` | One triplet repeated `L−2` times: `(L−2)(L−3)/2 / (L−2)` [2]. |
| INV-04 | Higher `S` ⇒ lower complexity | Repeated triplets increase `Σ c(c−1)/2` [1][2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` / `string` | required | Sequence to score | Upper-cased for the `string` overload; null `DnaSequence` throws |
| wordSize | `int` | 3 | Word (k-mer) size | ≥ 1; only k = 3 is source-defined [1][3] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| score | `double` | DUST score `Σ c(c−1)/2 / (L − wordSize + 1)`; ≥ 0; higher ⇒ lower complexity |

### 3.3 Preconditions and Validation

`string` input is upper-cased (T↔U not performed; DNA alphabet assumed). Null `DnaSequence` ⇒ `ArgumentNullException`; null/empty `string` ⇒ 0. `wordSize < 1` ⇒ `ArgumentOutOfRangeException`. A sequence shorter than one word (`L < wordSize`) ⇒ 0 (no word exists; defined-output convention, not a source value).

## 4. Algorithm

### 4.1 High-Level Steps

1. If the sequence is shorter than `wordSize`, return 0.
2. Slide a window of `wordSize` one base at a time and tally the count `c_t` of each distinct word over the `L − wordSize + 1` positions.
3. Sum `c_t·(c_t − 1)/2` over all distinct words.
4. Divide by the number of words (`L − wordSize + 1`, equal to `L − 2` for triplets) and return.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Word size:** k = 3 (triplets), hardcoded in DUST/SDUST [1][3]; exposed here as a parameter defaulting to 3.
- **Mask threshold:** 2.0 (reference default level `T = 20`, with `rw·10 > L·T` ⇔ score > 2.0) [3]; used by `MaskLowComplexity`.
- **Window size:** 64 bases (default for windowed masking) [1][2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateDustScore` | O(L·wordSize) | O(min(L, 4^wordSize)) | One pass; substring of length wordSize per position; hash map of distinct words |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateDustScore(DnaSequence, int)`: scores a `DnaSequence`.
- `SequenceComplexity.CalculateDustScore(string, int)`: scores a raw string (upper-cased).
- `SequenceComplexity.MaskLowComplexity(...)`: masks 64-base windows whose triplet score exceeds the 2.0 threshold (uses the same score core).

### 5.2 Current Behavior

The score is computed exactly as the formula in §2.2: numerator `Σ c(c−1)/2`, divisor = number of words (`L − wordSize + 1`). No suffix tree is used: the operation is a single linear scan tallying overlapping words, not an occurrence-enumeration query against a fixed text, so the repository suffix tree does not fit (recorded per the reuse policy). The earlier code divided by `(words − 1)`, which over-scaled the result; that was corrected to divide by the word count to conform to the `1/(L−2)` normalization [2][3].

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Numerator `Σ_t c_t·(c_t − 1)/2` over overlapping triplets [2][3].
- Normalization by the number of triplets `L − 2` (generalized to `L − wordSize + 1`) [2].
- Default triplet word size and mask threshold 2.0 (level 20) [1][3].

**Intentionally simplified:**

- General `wordSize` parameter: the score generalizes to non-triplet words by dividing by `L − wordSize + 1`; **consequence:** only `wordSize = 3` matches the published DUST/SDUST definition, other values are an extrapolation.

**Not implemented:**

- The symmetric perfect-interval masking rule (SDUST) that selects maximal high-scoring subwindows; **users should rely on:** `MaskLowComplexity`, which applies a simpler fixed-window threshold scan rather than the SDUST interval optimization.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `L < wordSize` ⇒ 0 | Assumption | No score defined by sources for sub-word input | accepted | Defined-output convention |
| 2 | General `wordSize` | Assumption | Only k=3 source-backed | accepted | Default 3; see 5.3 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `DnaSequence` | `ArgumentNullException` | Validation contract |
| Null/empty `string` | 0 | No words ⇒ minimal complexity |
| `L < wordSize` | 0 | No word exists (convention) |
| All-distinct triplets | 0 | INV-02 |
| Homopolymer length L | `(L−3)/2` | INV-03 |

### 6.2 Limitations

Heuristic, not a probabilistic significance test. Only triplet scoring (k = 3) is source-defined. The mask routine uses a fixed-window threshold scan rather than the symmetric SDUST interval rule, so masked boundaries may differ from `dustmasker`/`sdust`.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** `AAAAAA` (L = 6). Triplets at positions 0–3 are all `AAA` ⇒ `c_AAA = 4`. Numerator = `4·3/2 = 6`. Divisor = `L − 2 = 4`. Score = `6 / 4 = 1.5`.

`ACGTACGT` (L = 8): `ACG=2, CGT=2, GTA=1, TAC=1` ⇒ numerator = `1 + 1 = 2`, divisor = 6, score = `1/3 ≈ 0.3333`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceComplexity_CalculateDustScore_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateDustScore_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [SEQ-COMPLEX-DUST-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md)
- Related algorithms: [K-mer_Entropy](./K-mer_Entropy.md)

## 8. References

1. Morgulis A, Gertz EM, Schäffer AA, Agarwala R. 2006. A fast and symmetric DUST implementation to mask low-complexity DNA sequences. Journal of Computational Biology 13(5):1028–1040. https://doi.org/10.1089/cmb.2006.13.1028
2. Li H. 2025. Finding low-complexity DNA sequences with longdust. arXiv:2509.07357. https://arxiv.org/pdf/2509.07357
3. Li H. sdust — Symmetric DUST reference C implementation. https://raw.githubusercontent.com/lh3/sdust/master/sdust.c
