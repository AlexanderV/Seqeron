# DUST Score

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-DUST-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Complete (score + SDUST masking) |
| Last Reviewed | 2026-09-28 |

## 1. Overview

The DUST score is a heuristic measure of nucleotide-sequence low complexity introduced for masking repetitive DNA in BLAST [1]. It counts how often each overlapping triplet (3-mer) recurs in a sequence and combines those counts into a single score: highly repetitive sequences (e.g. homopolymers, simple satellites) score high, while sequences with all-distinct triplets score 0. It is a deterministic heuristic, not a probabilistic model; this unit implements the per-sequence/per-window complexity score, the quantity DUST/SDUST thresholds to decide masking.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Low-complexity regions (homopolymers, tandem repeats, simple sequence repeats) produce spurious alignments and inflate database-search noise. DUST ("Decreasing Uncertainty in Sequence Tags") masks such regions by scoring fixed-size windows and masking windows whose score exceeds a threshold [1]. The score function is unchanged between the original and the symmetric (SDUST) implementation; only the masking rule differs [1].

### 2.2 Core Model

For a sequence `x` of length `L`, let `c_t` be the number of occurrences of triplet `t` among the `ℓ = L − 2` overlapping triplets. The complexity score is [1][3][4]:

```
S(x) = ( Σ_t  c_t·(c_t − 1)/2 ) / (ℓ − 1)
```

The numerator `Σ_t c_t(c_t−1)/2` counts pairs of identical triplets (the reference accumulates it as `r += c[t]++` [3]). The normaliser is `ℓ − 1`: NCBI dustmasker tests `10·r > thresholds_[ℓ−1]` with `thresholds_[i] = i·level` [4], and lh3/sdust tests `new_r·10 > T·new_l` with `new_l = kdq_size(w) − i − 1 = ℓ − 1` [3]. Numerically confirmed with the compiled sdust binary (`-t 20`): an isolated 7-A run (ℓ = 5, r = 10, 10/4 = 2.5) is masked while a 6-A run (6/3 = 2.0) is not; `(AC)×6` (20/9) is masked while `ACACACACACA` (16/8 = 2.0) is not. A divisor of `ℓ` would mask none of them. (The longdust README/preprint [2] writes the SDUST score with `ℓ(x)`; this does not match the reference code and is not used.)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `S(x) ≥ 0` | Each `c(c−1)/2 ≥ 0` and the divisor `ℓ − 1 > 0` for `L ≥ 4`; `ℓ ≤ 1` ⇒ 0. |
| INV-02 | All-distinct triplets ⇒ `S(x) = 0` | Every `c_t = 1` ⇒ `c(c−1)/2 = 0` [2]. |
| INV-03 | Homopolymer of length `L ≥ 4` ⇒ `S = (L−2)/2` | One triplet repeated `ℓ = L−2` times: `(L−2)(L−3)/2 / (L−3)`. |
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
| score | `double` | DUST score `Σ c(c−1)/2 / (ℓ − 1)`, `ℓ = L − wordSize + 1`; ≥ 0; higher ⇒ lower complexity |

### 3.3 Preconditions and Validation

`string` input is upper-cased (T↔U not performed; DNA alphabet assumed). Null `DnaSequence` ⇒ `ArgumentNullException`; null/empty `string` ⇒ 0. `wordSize < 1` ⇒ `ArgumentOutOfRangeException`. Fewer than two words (`ℓ ≤ 1`) ⇒ 0 (no word pair exists; defined-output convention, not a source value). `MaskLowComplexity`: `windowSize < 3` or a negative/NaN/infinite `threshold` ⇒ `ArgumentOutOfRangeException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. If fewer than two words exist (`ℓ = L − wordSize + 1 ≤ 1`), return 0.
2. Tally the count `c_t` of each overlapping word with the canonical `KmerAnalyzer.CountKmers`.
3. Sum `c_t·(c_t − 1)/2` over all distinct words.
4. Divide by `ℓ − 1` and return.

**SDUST masking (`MaskLowComplexity`)** — a port of lh3/sdust `sdust_core` [3]: slide a window of `W` bases (≤ W − 2 triplets); maintain the window score `rw` and the score `rv` of the longest window suffix in which no triplet occurs more than `2·threshold` times; when `rw > threshold·L_suffix`, extend leftwards from that suffix to find every *perfect interval* (score > threshold and ≥ the best score of any perfect interval it contains); intervals leaving the window are emitted and merged when overlapping/adjacent. Non-ACGT characters split the scan.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Word size:** k = 3 (triplets), hardcoded in DUST/SDUST [1][3]; exposed here as a parameter defaulting to 3.
- **Mask threshold:** 2.0 (reference default level `T = 20`, with `rw·10 > L·T` ⇔ score > 2.0) [3]; used by `MaskLowComplexity`.
- **Window size:** 64 bases (default for windowed masking) [1][2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateDustScore` | O(L·wordSize) | O(min(L, 4^wordSize)) | One pass via `KmerAnalyzer.CountKmers` |
| `MaskLowComplexity` | O(L·W) typical, O(L·W³) worst | O(W + 64) | As lh3/sdust |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateDustScore(DnaSequence, int)`: scores a `DnaSequence`.
- `SequenceComplexity.CalculateDustScore(string, int)`: scores a raw string (upper-cased).
- `SequenceComplexity.MaskLowComplexity(...)`: symmetric DUST (SDUST) masking of perfect intervals; default W = 64, threshold 2.0 (level 20).

### 5.2 Current Behavior

The score is computed exactly as the formula in §2.2: numerator `Σ c(c−1)/2`, divisor `ℓ − 1`. Word counting delegates to `KmerAnalyzer.CountKmers`. History: the original code divided by `words − 1` (correct); a 2026-06 change "corrected" it to the word count from a misreading of the longdust restatement [2]; the 2026-09 review restored `ℓ − 1` after confirming it against the NCBI symdust and lh3/sdust sources and the compiled sdust binary. `MaskLowComplexity` was previously a fixed 64-bp window scan that masked whole windows (and skipped sequences shorter than the window); it is now the real SDUST algorithm, verified identical to the lh3/sdust binary on 3,000 random repeat-rich sequences (W ∈ {3,5,8,16,30,64,100}, T ∈ {10,12,15,20,25,30}) and a 1-Mb sequence.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Numerator `Σ_t c_t·(c_t − 1)/2` over overlapping triplets [1][3].
- Normalization by `ℓ − 1` [1][3][4] (generalized to `L − wordSize`).
- Default triplet word size, window 64 and mask threshold 2.0 (level 20) [3][4].
- Symmetric perfect-interval masking (SDUST) [1][3].

**Intentionally simplified:**

- General `wordSize` parameter of `CalculateDustScore`: generalizes to non-triplet words by dividing by `ℓ − 1 = L − wordSize`; **consequence:** only `wordSize = 3` matches the published DUST/SDUST definition, other values are an extrapolation. `MaskLowComplexity` always uses triplets.

**Not implemented:** dustmasker's `linker` post-merge of intervals separated by ≤ 1 bp (sdust has none either; output equals sdust).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `ℓ ≤ 1` ⇒ 0 | Assumption | No score defined by sources (ℓ − 1 = 0) | accepted | Defined-output convention |
| 2 | General `wordSize` | Assumption | Only k=3 source-backed | accepted | Default 3; see 5.3 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `DnaSequence` | `ArgumentNullException` | Validation contract |
| Null/empty `string` | 0 | No words ⇒ minimal complexity |
| `ℓ ≤ 1` (L ≤ 3 for triplets) | 0 | No word pair exists (convention) |
| All-distinct triplets | 0 | INV-02 |
| Homopolymer length L ≥ 4 | `(L−2)/2` | INV-03 |

### 6.2 Limitations

Heuristic, not a probabilistic significance test. Only triplet scoring (k = 3) is source-defined. Masked intervals equal lh3/sdust output; dustmasker differs from sdust only in assembly gaps (sdust README) — not reachable through `DnaSequence`, which admits only A/C/G/T.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** `AAAAAA` (L = 6). Triplets at positions 0–3 are all `AAA` ⇒ `c_AAA = 4`, ℓ = 4. Numerator = `4·3/2 = 6`. Divisor = `ℓ − 1 = 3`. Score = `6 / 3 = 2.0`.

`ACGTACGT` (L = 8, ℓ = 6): `ACG=2, CGT=2, GTA=1, TAC=1` ⇒ numerator = `1 + 1 = 2`, divisor = 5, score = `0.4`.

Masking: `ACGTTGCAGTCATGCGATC` + `A×7` + `TGCATCGGATCCTAGGCTAA` ⇒ sdust/`MaskLowComplexity` mask [19, 26).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceComplexity_CalculateDustScore_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceComplexity_CalculateDustScore_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [SEQ-COMPLEX-DUST-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md)
- Related algorithms: [K-mer_Entropy](./K-mer_Entropy.md)

## 8. References

1. Morgulis A, Gertz EM, Schäffer AA, Agarwala R. 2006. A fast and symmetric DUST implementation to mask low-complexity DNA sequences. Journal of Computational Biology 13(5):1028–1040. https://doi.org/10.1089/cmb.2006.13.1028
2. Li H. 2025. Finding low-complexity DNA sequences with longdust. arXiv:2509.07357. https://arxiv.org/pdf/2509.07357
3. Li H. sdust — Symmetric DUST reference C implementation. https://raw.githubusercontent.com/lh3/sdust/master/sdust.c
4. NCBI C++ Toolkit — dustmasker symmetric DUST (`src/algo/dustmask/symdust.cpp`, `include/algo/dustmask/symdust.hpp`). https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/dustmask/symdust.cpp
