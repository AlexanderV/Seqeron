# Lempel–Ziv Complexity

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-COMPRESS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

The Lempel–Ziv (1976) complexity measures how compressible a finite sequence is by counting the number of distinct components produced when the sequence is parsed left-to-right, each new component being the shortest substring not yet encountered in the scanned history [1][2]. It is a deterministic, combinatorial complexity measure (no probabilistic model): repetitive sequences yield few components (low complexity), while diverse/random sequences yield many [1]. The repository exposes a raw component count, a length-normalized variant, and a `EstimateCompressionRatio` entry point that returns the normalized value, used to flag low-complexity / repetitive regions in nucleotide sequences.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A finite sequence over a finite alphabet is "built up" by two operations: copying a symbol/substring already seen, or generating a genuinely new substring. The complexity is the number of new substrings (components) needed to reconstruct the sequence by this exhaustive-history production process [1][2]. The measure underlies the LZ77/LZ78/LZW compressors and is widely reused in bio-sequence analysis to quantify repetitiveness [5].

### 2.2 Core Model

Parse the sequence `S` of length `n` left-to-right. Maintain the set of components already produced. Starting at index `ind` with length `inc = 1`, let `sub = S[ind .. ind+inc)`. If `sub` is already a produced component, extend it (`inc ← inc + 1`); otherwise emit `sub` as a new component, advance `ind ← ind + inc`, and reset `inc ← 1`. The Lempel–Ziv complexity `c(S)` is the number of components produced [2][3].

Normalization removes the length dependence of `c` [4][5]: with `b` the alphabet size (number of distinct symbols present) and `b(n) = n / log_b(n)` the asymptotic upper bound for a uniformly random sequence [6],

`LZ_norm = c / (n / log_b(n))`.

`LZ_norm → 1` for a maximally complex (random) sequence; smaller values indicate more compressible input [6].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `c(S) = 0` iff `S` is empty | no component is produced for an empty scan [3] |
| INV-02 | `c(S) ≥ 1` for any non-empty `S` | the first symbol is always a new component [3] |
| INV-03 | `c(S) ≤ n` for `|S| = n` | each component has length ≥ 1 [1][3] |
| INV-04 | a homopolymer has strictly lower `c` than an all-distinct string of equal length | productivity buildup: a run reuses components, distinct symbols each start a new one [1][2] |
| INV-05 | `EstimateCompressionRatio(S) = LZ_norm(S)` | implemented as a thin delegate (design) |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Lempel–Ziv complexity | Shannon entropy |
|--------|-----------------------|-----------------|
| Captures | sequential / positional repetition | symbol-frequency distribution only |
| Sensitive to order | yes | no |
| Model | combinatorial, model-free [1] | probabilistic (per-symbol frequencies) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` / `string` | required | sequence to analyze | upper-cased internally; any alphabet accepted |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CalculateLempelZivComplexity` | `int` | raw component count `c(S)` |
| `CalculateNormalizedLempelZivComplexity` | `double` | `c / (n / log_b(n))`; raw count if `b < 2` |
| `EstimateCompressionRatio` | `double` | normalized Lempel–Ziv complexity (delegates to the above) |

### 3.3 Preconditions and Validation

Null `DnaSequence` throws `ArgumentNullException`. Null/empty `string` returns 0 (raw) or 0 (normalized). Input is upper-cased (`ToUpperInvariant`); parsing is alphabet-agnostic (works for DNA `{A,C,G,T}` and the binary examples). For normalization, `b` is the number of distinct symbols actually present; when `b < 2` (single distinct symbol, so `log_b n` is undefined) the raw count is returned [4].

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case the input; empty → 0.
2. Exhaustive-history parse: grow the running substring while it is a known component; otherwise emit it as a new component and restart at length 1.
3. Return the number of distinct components (raw complexity).
4. For normalization, compute `b` (distinct symbols), `b(n) = n / log_b(n)`, and return `c / b(n)`; if `b < 2` return `c`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| raw complexity | O(n) amortized over the scan; substring/hash work makes it O(n²) worst case for the hash-set variant | O(n) | one left-to-right pass with a component set [2][3] |
| normalization | O(n) | O(σ) | adds a distinct-symbol count |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateLempelZivComplexity(DnaSequence|string)`: raw LZ76 component count.
- `SequenceComplexity.CalculateNormalizedLempelZivComplexity(DnaSequence|string)`: length-normalized complexity.
- `SequenceComplexity.EstimateCompressionRatio(DnaSequence|string)`: delegates to the normalized method (registry-canonical name).

### 5.2 Current Behavior

The parse uses a `HashSet<string>` of produced components and the set-based exhaustive-history rule from the Naereen reference implementation [3]: a trailing partial substring that never becomes a new component on its own is NOT counted (it was never added to the set). This is the convention whose worked values are reproducible and is the contract enforced by the tests. Suffix-tree reuse was evaluated and **not** used: LZ76 parsing is a single left-to-right factorization, not a multi-query exact-match search, so the repository suffix tree does not fit; a single linear scan with a component set is the correct structure.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exhaustive-history left-to-right component parse and raw count `c(S)` [1][2][3].
- Length normalization `c / (n / log_b(n))` with `b` = alphabet size [4][5][6].

**Intentionally simplified:**

- (none)

**Not implemented:**

- LZ77/LZ78 dictionary *encoding/decoding* (the compressed bitstream); only the complexity *count* is computed. Users needing an actual compressor should rely on a dedicated compression library.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | trailing-partial-component not counted | Assumption | last component count may differ by ±1 vs the Wikipedia pseudocode `if v≠1` rule | accepted | set-based reference convention [3]; worked values reproducible |
| 2 | normalization log base = distinct symbols present; `b<2` ⇒ raw count | Assumption | normalized value undefined for single-symbol input | accepted | [4]; documented degenerate handling |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| empty / null string | 0 | INV-01 [3] |
| null `DnaSequence` | `ArgumentNullException` | sibling-method convention |
| single base `"A"` | raw 1 | INV-02 [3] |
| homopolymer `"0"×16` | raw 5 (`0/00/000/0000/00000`) | exhaustive-history parse [3] |
| single-symbol input (normalized) | raw count (log base undefined) | `b < 2` rule [4] |

### 6.2 Limitations

Raw complexity grows with sequence length, so only the normalized value is comparable across sequences of different lengths [4][5]. The measure is not an actual compressed size; it does not model nucleotide-specific biology. The hash-set parse is O(n²) worst case on highly repetitive long inputs.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
int c = SequenceComplexity.CalculateLempelZivComplexity("1001111011000010"); // 8
double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("1001111011000010"); // 2.0
```

**Numerical walk-through:** `1001111011000010` parses as `1 / 0 / 01 / 11 / 10 / 110 / 00 / 010` → 8 components [3]. Normalized: `n=16`, `b=2`, `log₂16 = 4`, `b(n) = 16/4 = 4`, `LZ_norm = 8/4 = 2.0` [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceComplexity_EstimateCompressionRatio_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_EstimateCompressionRatio_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [SEQ-COMPLEX-COMPRESS-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md)

## 8. References

1. Lempel, A., & Ziv, J. 1976. On the Complexity of Finite Sequences. IEEE Transactions on Information Theory 22(1):75–81. https://doi.org/10.1109/TIT.1976.1055501
2. Wikipedia. Lempel–Ziv complexity. https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity
3. Besson, L. (Naereen). Lempel-Ziv_Complexity (Python reference implementation). https://github.com/Naereen/Lempel-Ziv_Complexity/blob/master/src/lempel_ziv_complexity.py
4. Vallat, R. entropy / AntroPy — `lziv_complexity`. https://raphaelvallat.com/entropy/build/html/generated/entropy.lziv_complexity.html
5. Zhang, Y., Hao, J., Zhou, C., & Chang, K. 2009. Normalized Lempel-Ziv complexity and its application in bio-sequence analysis. Journal of Mathematical Chemistry 46(4):1203–1212. https://doi.org/10.1007/s10910-008-9512-2
6. Hu, J., Gao, J., & Principe, J.C. 2006. Analysis of biomedical signals by the Lempel-Ziv complexity. arXiv:nlin/0608049. https://arxiv.org/abs/nlin/0608049
