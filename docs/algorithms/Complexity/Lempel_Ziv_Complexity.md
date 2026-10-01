# Lempel–Ziv Complexity

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-COMPRESS-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

The Lempel–Ziv (1976) complexity measures how compressible a finite sequence is by counting the components of its *exhaustive history*: parsing left-to-right, each component is extended while it can still be copied from the text before it (overlap allowed) and closed by the first symbol that makes it new [1][2]. It is a deterministic, combinatorial complexity measure (no probabilistic model): repetitive sequences yield few components (low complexity), while diverse/random sequences yield many [1]. The repository exposes a raw component count, a length-normalized variant, and a `EstimateCompressionRatio` entry point that returns the normalized value, used to flag low-complexity / repetitive regions in nucleotide sequences.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A finite sequence over a finite alphabet is "built up" by two operations: copying a symbol/substring already seen, or generating a genuinely new substring. The complexity is the number of new substrings (components) needed to reconstruct the sequence by this exhaustive-history production process [1][2]. The measure underlies the LZ77/LZ78/LZW compressors and is widely reused in bio-sequence analysis to quantify repetitiveness [5].

### 2.2 Core Model

**LZ76 exhaustive history** [1][6]. Parse `S` (length `n`) left-to-right. A component starting at position `p` is `w = S[p .. p+L)` with `L` the smallest length such that `w` is **not** a substring of `S[0 .. p+L−1)` (i.e. `w` without its last symbol is reproducible by copying from some start `< p`, overlap allowed, but `w` itself is not); if the end of `S` is reached while `w` is still reproducible, the remainder is the last component. `c(S)` is the number of components. Examples: `0001101001000101 = 0·001·10·100·1000·101` (c = 6, Lempel & Ziv 1976 [1]); `1001111011000010 = 1/0/01/1110/1100/0010` (c = 6 [2][4]); `0×16 = 0/0…0` (c = 2).

The reference algorithm is the scan of Kaspar & Schuster (1987) [3] (antropy `_lz_complexity` [4]). Production computes the same factorization from the Longest-Previous-Factor array `LPF[q] = max_{j<q} lcp(S[q..], S[j..])` [8] (component at `q` = `S[q .. q+LPF[q]]`), which is value-identical to the scan (checked against antropy and a brute-force implementation of the definition on 20 000 random strings) but runs in O(n log² n).

**Not LZ78.** The "set of previously seen phrases" parse (each new phrase = a seen phrase + one symbol; e.g. Naereen `lempel_ziv_complexity` [7]) is the Ziv–Lempel 1978 incremental parsing, not the LZ76 measure: it gives 8 for `1001111011000010` and 5 for `0×16`. It was used here before 2026-09 and was replaced.

Normalization removes the length dependence of `c` [4][5]: with `b` the alphabet size (number of distinct symbols present) and `b(n) = n / log_b(n)` the asymptotic upper bound for a uniformly random sequence [6],

`LZ_norm = c / (n / log_b(n))`.

`LZ_norm → 1` for a maximally complex (random) sequence; smaller values indicate more compressible input [6].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `c(S) = 0` iff `S` is empty | no component is produced for an empty scan |
| INV-02 | `c(S) ≥ 1` for any non-empty `S` | the first symbol is always a new component [1] |
| INV-03 | `c(S) ≤ n` for `|S| = n` | each component has length ≥ 1 [1] |
| INV-04 | a homopolymer has `c = min(n, 2)`, the minimum for its length; strictly lower than an all-distinct string of equal length n ≥ 3 | the whole run after the first symbol is one self-overlapping copy [1] |
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
| `CalculateNormalizedLempelZivComplexity` | `double` | `c / (n / log_b(n))`, `b := max(b, 2)`; raw count (1) if `n = 1` |
| `EstimateCompressionRatio` | `double` | normalized Lempel–Ziv complexity (delegates to the above) |

### 3.3 Preconditions and Validation

Null `DnaSequence` throws `ArgumentNullException`. Null/empty `string` returns 0 (raw) or 0 (normalized). Input is upper-cased (`ToUpperInvariant`); parsing is alphabet-agnostic (works for DNA `{A,C,G,T}` and the binary examples). For normalization, `b` is the number of distinct symbols actually present, clamped to `b := max(b, 2)` as in antropy [4]; for `n = 1` (`log_b 1 = 0`, where antropy raises a division-by-zero) the raw count 1 is returned.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case the input; empty → 0.
2. Build the suffix array (prefix doubling), the Kasai LCP array, and from them the Longest-Previous-Factor array `LPF` (nearest smaller text position on either side in suffix-array order) [8]. Starting at `q = 0`: count a component and jump `q ← q + LPF[q] + 1` until `q ≥ n` (a final still-reproducible remainder counts as one component).
3. Return the component count (raw complexity).
4. For normalization, compute `b = max(distinct symbols, 2)`, `b(n) = n / log_b(n)`, and return `c / b(n)` (`n = 1` → `c`).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| raw complexity | O(n log² n) | O(n) | suffix array + LCP + LPF [8]; the Kaspar–Schuster scan [3] would be O(n²/log n) on random input |
| normalization | O(n) | O(σ) | adds a distinct-symbol count |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateLempelZivComplexity(DnaSequence|string)`: raw LZ76 (exhaustive-history) component count.
- `SequenceComplexity.CalculateNormalizedLempelZivComplexity(DnaSequence|string)`: length-normalized complexity.
- `SequenceComplexity.EstimateCompressionRatio(DnaSequence|string)`: delegates to the normalized method (registry-canonical name).

### 5.2 Current Behavior

The raw count is computed from the LPF array [8]; it reproduces the Kaspar–Schuster scan [3] (antropy `_lz_complexity` [4], Wikipedia pseudocode [2]) exactly, including the trailing reproducible component (`if len ≠ 1 then c += 1`). Values are locked against antropy doctests, the Lempel–Ziv (1976) and Estévez-Rams worked examples, and an antropy-computed random-DNA dataset; the fuzz suite cross-checks against a brute-force implementation of the definition.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- LZ76 exhaustive-history component count `c(S)` [1], computed from the LPF array [8] — value-identical to the Kaspar–Schuster scan [3].
- Length normalization `c / (n / log_b(n))` with `b` = alphabet size [4][5][6].

**Intentionally simplified:**

- (none)

**Not implemented:**

- LZ77/LZ78 dictionary *encoding/decoding* (the compressed bitstream); only the complexity *count* is computed. Users needing an actual compressor should rely on a dedicated compression library.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | normalization log base = distinct symbols present, clamped to ≥ 2 | Convention | matches antropy [4] | accepted | Zhang et al. 2009 [5] |
| 2 | `n = 1` ⇒ normalized returns raw count 1 | Convention | formula undefined (`log_b 1 = 0`; antropy raises) | accepted | documented degenerate handling |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| empty / null string | 0 | INV-01 |
| null `DnaSequence` | `ArgumentNullException` | sibling-method convention |
| single base `"A"` | raw 1, normalized 1 | INV-02; `n = 1` guard |
| homopolymer `"0"×16` | raw 2 (`0/0…0`), normalized `2/(16/log₂16) = 0.5` | exhaustive history [1]; antropy [4] |
| single-symbol input (normalized) | base clamped to 2 | antropy `base = 2 if base < 2` [4] |

### 6.2 Limitations

Raw complexity grows with sequence length, so only the normalized value is comparable across sequences of different lengths [4][5]. The measure is not an actual compressed size; it does not model nucleotide-specific biology. Time is O(n log² n) (200 kb random DNA in ≈ 1 s).

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
int c = SequenceComplexity.CalculateLempelZivComplexity("1001111011000010"); // 6
double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity("1001111011000010"); // 1.5
```

**Numerical walk-through:** `1001111011000010` has exhaustive history `1 / 0 / 01 / 1110 / 1100 / 0010` → 6 components [2][4] (e.g. `1110`: `111` is copyable from earlier text with overlap, `1110` is not). Normalized: `n=16`, `b=2`, `log₂16 = 4`, `b(n) = 16/4 = 4`, `LZ_norm = 6/4 = 1.5` (antropy doctest [4]).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceComplexity_EstimateCompressionRatio_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceComplexity_EstimateCompressionRatio_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [SEQ-COMPLEX-COMPRESS-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md)

## 8. References

1. Lempel, A., & Ziv, J. 1976. On the Complexity of Finite Sequences. IEEE Transactions on Information Theory 22(1):75–81. https://doi.org/10.1109/TIT.1976.1055501
2. Wikipedia. Lempel–Ziv complexity. https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity
3. Kaspar, F., & Schuster, H. G. 1987. Easily calculable measure for the complexity of spatiotemporal patterns. Physical Review A 36(2):842–848. https://doi.org/10.1103/PhysRevA.36.842
4. Vallat, R. AntroPy 0.2.2 — `lziv_complexity` / `_lz_complexity` (source `antropy/entropy.py`, PyPI wheel). https://github.com/raphaelvallat/antropy
5. Zhang, Y., Hao, J., Zhou, C., & Chang, K. 2009. Normalized Lempel-Ziv complexity and its application in bio-sequence analysis. Journal of Mathematical Chemistry 46(4):1203–1212. https://doi.org/10.1007/s10910-008-9512-2
6. Hu, J., Gao, J., & Principe, J.C. 2006. Analysis of biomedical signals by the Lempel-Ziv complexity. arXiv:nlin/0608049. https://arxiv.org/abs/nlin/0608049 ; Estévez-Rams, E. et al. 2013. On the non-randomness of maximum Lempel Ziv complexity sequences of finite size. arXiv:1311.0546 (exhaustive-history definition; example `010011101101100 = 0.1.00.11.101.101100`, C = 6).
7. Besson, L. (Naereen). Lempel-Ziv_Complexity 0.2.2 (PyPI) — implements the LZ78-style incremental parse, NOT LZ76 (kept as a counter-reference).
8. Crochemore, M., & Ilie, L. 2008. Computing Longest Previous Factor in linear time and applications. Information Processing Letters 106(2):75–80. https://doi.org/10.1016/j.ipl.2007.10.006
