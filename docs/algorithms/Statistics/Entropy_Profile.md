# Shannon Entropy Profile

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-ENTROPY-PROFILE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

The Shannon entropy profile measures local sequence complexity by computing the Shannon entropy H = −Σ pᵢ log₂ pᵢ (in bits) over each fixed-width window slid along the sequence [1][2]. Each window contributes one value; low values mark low-complexity / repetitive regions (a homopolymer window is 0 bits) and high values mark compositionally diverse regions (a uniform 4-nucleotide window is the 2-bit maximum) [2][3]. The computation is exact and deterministic: it is the closed-form information-theoretic entropy of the per-symbol frequency distribution of each window. It is used to scan genomes for low-complexity tracts and compositional shifts [3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Shannon entropy quantifies the average uncertainty (choice/information) of a discrete source [1]. Applied to a biological sequence window, the "source" is the empirical distribution of symbols in that window; entropy then measures how evenly the symbols are distributed, i.e. local compositional complexity [3].

### 2.2 Core Model

For a window with symbol frequencies pᵢ (counts divided by the number of counted symbols), Shannon entropy is

H = −Σᵢ pᵢ log_b pᵢ [1][2].

With base b = 2 the unit is bits (shannons) [2]. Terms with pᵢ = 0 contribute 0 (0·log 0 ≡ 0) [1][2]. For k equally-likely symbols H attains its maximum log₂ k; for the 4-letter DNA alphabet the maximum is log₂ 4 = 2 bits [2][3]. The profile is the ordered sequence of H values, one per window of width W slid along the sequence [3] (IntechOpen Eq. 3: yᵢ = −Σⱼ pᵢⱼ log pᵢⱼ).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each profile value H ≥ 0 | every term −pᵢ log₂ pᵢ ≥ 0 for pᵢ ∈ (0,1] [2] |
| INV-02 | Each profile value H ≤ log₂ k (k = distinct symbols in the window; ≤ 2 bits for DNA) | maximum entropy at the uniform distribution is log₂ k [2][3] |
| INV-03 | A homopolymer window yields H = 0 | single symbol has p = 1, log₂ 1 = 0 [1][2] |
| INV-04 | A window with all symbols equally frequent yields H = log₂ k | uniform distribution attains the maximum [2] |
| INV-05 | Number of windows = ⌊(n − W)/step⌋ + 1 when W ≤ n, else 0 | windows start at offsets 0, step, 2·step, … ≤ n − W [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | sequence to profile | letters are counted (case-folded); non-letters ignored |
| `windowSize` | `int` | 50 | window width W in symbols | a window is produced only when W ≤ length |
| `stepSize` | `int` | 1 | window advance in symbols | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IEnumerable<double>` | one Shannon entropy value (bits) per window, in offset order |

### 3.3 Preconditions and Validation

Null or empty `sequence`, or `windowSize` greater than the sequence length, yields an empty profile (no exception). Counting is case-insensitive (input is upper-cased) and restricted to letters; degenerate/`N` symbols are counted as their own symbol. There is no T↔U normalization — U and T are distinct symbols if both appear. Indexing of window offsets is 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. If `sequence` is null/empty or `windowSize` > length, produce nothing.
2. For each offset i = 0, stepSize, 2·stepSize, … while i ≤ length − windowSize:
3. Take the window `sequence[i .. i+windowSize)`, count its letter frequencies (case-folded), and compute H = −Σ pᵢ log₂ pᵢ.
4. Yield H for that window.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Full profile | O(n × W) | O(σ) | n = length, W = window width, σ = alphabet size; one O(W) pass per window |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateEntropyProfile(string, int windowSize=50, int stepSize=1)`: sliding-window driver; yields one entropy per window.
- `SequenceStatistics.CalculateShannonEntropy(string)`: per-window kernel; computes H = −Σ pᵢ log₂ pᵢ over the letter frequencies of the argument.

### 5.2 Current Behavior

The profile is a lazy `IEnumerable<double>` (deferred, streaming). Each window is materialized as a substring and delegated to `CalculateShannonEntropy`, which case-folds, counts only `char.IsLetter` symbols, and uses `Math.Log2` (base 2 → bits). No suffix tree is used: this is not a substring-search/occurrence problem but a per-window frequency-counting computation, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- H = −Σ pᵢ log₂ pᵢ in bits, base-2 logarithm [1][2].
- Zero-probability convention (absent symbols contribute 0) [1][2].
- Maximum of log₂ k at the uniform distribution; 2-bit ceiling for the DNA alphabet [2][3].
- Per-window sliding computation of the entropy [3].

**Intentionally simplified:**

- Mono-symbol (k=1) alphabet for pᵢ: entropy is taken over single-symbol frequencies, not k-mer / block entropy. **Consequence:** the profile reflects single-symbol composition only; higher-order correlations between adjacent symbols are not captured.

**Not implemented:**

- Higher-order (block / n-mer) entropy; **users should rely on:** k-mer based complexity measures such as `SequenceComplexity.CalculateKmerEntropy` for correlation-aware analysis.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | empty profile | guarded input (§3.3) |
| windowSize > length | empty profile | no full window exists (INV-05) |
| windowSize == length | single value | exactly one window |
| homopolymer window | 0.0 bits | INV-03 [1][2] |
| uniform DNA window | 2.0 bits | INV-04/INV-02 [2][3] |

### 6.2 Limitations

Single-symbol entropy does not distinguish sequences with identical composition but different order (e.g., `AATT` and `ATAT` both yield 1 bit). It is alphabet-sensitive: protein windows can exceed 2 bits since k > 4. No statistical significance / background-model correction is applied.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Sliding Shannon entropy (bits), window 4, step 1, over "AAATGC".
var profile = SequenceStatistics.CalculateEntropyProfile("AAATGC", windowSize: 4, stepSize: 1).ToArray();
// profile == [0.8112781244591328, 1.5, 2.0]
//   AAAT -> -(3/4 log2 3/4 + 1/4 log2 1/4) = 0.8112781244591328
//   AATG -> -(1/2 log2 1/2 + 2*(1/4 log2 1/4)) = 1.5
//   ATGC -> log2 4 = 2.0
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateEntropyProfile_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateEntropyProfile_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [SEQ-ENTROPY-PROFILE-001-Evidence.md](../../../docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md)
- Related algorithms: [Dinucleotide_Analysis](./Dinucleotide_Analysis.md)

## 8. References

1. Shannon, C. E. 1948. A Mathematical Theory of Communication. Bell System Technical Journal, 27(3):379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x
2. Wikipedia contributors. Entropy (information theory). https://en.wikipedia.org/wiki/Entropy_(information_theory) (accessed 2026-06-14).
3. Entropy-Based Biological Sequence Study. IntechOpen. https://www.intechopen.com/chapters/75997 (accessed 2026-06-14).
