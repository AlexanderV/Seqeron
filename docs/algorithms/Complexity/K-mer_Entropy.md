# K-mer Entropy

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-KMER-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

K-mer entropy measures the sequence complexity of a DNA string as the Shannon entropy (in bits) of the frequency distribution of its overlapping k-mers. It quantifies how uniformly the k-mers are distributed: a low value indicates a few dominant k-mers (repeats, homopolymers, low complexity), while a high value indicates a near-uniform k-mer distribution (high complexity) [1]. The computation is exact and deterministic for a given sequence and k.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Detecting low-complexity DNA regions is a standard pre-processing step for sequence alignment and search. Modelling complexity through the entropy of k-mer (length-k substring) counts is the basis of tools such as `longdust` [1]; the entropy saturates at its maximum for random uniform sequences and drops toward zero for repetitive ones [2].

### 2.2 Core Model

Decompose a sequence of length L into its overlapping k-mers using a sliding window of step 1, giving N = L − k + 1 k-mers [1]. Let n_i be the count of the i-th distinct k-mer and p_i = n_i / N its relative frequency, so Σ p_i = 1. The Shannon entropy is

> H = − Σ_i p_i · log₂(p_i)   (bits) [1][3]

The base-2 logarithm yields entropy in bits [2]. This is the Shannon entropy H(X) = −Σ p(x) log p(x) of the k-mer distribution [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ H ≤ log₂(N), N = L − k + 1 | Shannon entropy is bounded below by 0 and above by log_b of the number of outcomes [3] |
| INV-02 | A single distinct k-mer (deterministic distribution) ⇒ H = 0 | H = 0 iff one outcome has p = 1 [3]; equivalently a fully repetitive sequence [1] |
| INV-03 | All k-mers distinct (uniform) ⇒ H = log₂(N) | Entropy is maximised, equal to log_b(n), under the uniform distribution [3] |
| INV-04 | Result is invariant to letter case | Input is normalised to upper-case (DnaSequence and the string overload) before counting |

### 2.5 Comparison with Related Methods

| Aspect | K-mer entropy (this) | Per-base Shannon entropy (`CalculateShannonEntropy`) |
|--------|----------------------|------------------------------------------------------|
| Alphabet of outcomes | distinct k-mers | the 4 nucleotides (k = 1 over A/C/G/T only) |
| Sensitive to local order | yes (k ≥ 2 captures di-/tri-nucleotide structure) | no (composition only) |
| Maximum value | log₂(L − k + 1) | 2 bits (log₂ 4) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` or `string` | required | Sequence to analyse | string is upper-cased; null/empty string → 0; null DnaSequence → throws |
| k | `int` | 2 | K-mer (window) length | k ≥ 1; k > L ⇒ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `double` | Shannon entropy in bits of the overlapping k-mer frequency distribution; 0 when L < k |

### 3.3 Preconditions and Validation

Indexing is 0-based over positions 0..L−k (inclusive). The accepted alphabet is unconstrained: every distinct length-k substring is treated as a symbol (no IUPAC filtering). Input is normalised to upper-case, so the result is case-insensitive. `k < 1` raises `ArgumentOutOfRangeException`; a null `DnaSequence` raises `ArgumentNullException`; a null/empty `string` returns 0; `k > L` returns 0 (no k-mers exist).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate k (≥ 1) and the sequence; normalise case.
2. If L < k, return 0.
3. Slide a window of length k by one position at a time, counting each distinct k-mer; total count N = L − k + 1.
4. For each distinct k-mer, p_i = n_i / N; accumulate −p_i · log₂(p_i).
5. Return the sum (bits).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateKmerEntropy | O(N · k) | O(D · k) | N = L − k + 1 windows; each substring build/hash is O(k); D = number of distinct k-mers stored in the dictionary |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateKmerEntropy(DnaSequence, int)`: canonical entry; validates and delegates to the core.
- `SequenceComplexity.CalculateKmerEntropy(string, int)`: string overload; upper-cases then delegates to the same core.
- `SequenceComplexity.CalculateKmerEntropyCore(string, int)` (private): counts overlapping k-mers in a dictionary and applies the entropy formula.

### 5.2 Current Behavior

K-mers are enumerated with a single linear scan and a `Dictionary<string,int>` of counts; entropy is then computed over the dictionary values. The repository suffix tree was evaluated and **not** used: this is a single linear pass building a full frequency table (every position visited once), not a repeated occurrence-query workload, so a suffix tree adds construction overhead without changing the linear cost or the required output.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Overlapping k-mer decomposition with N = L − k + 1 [1].
- H = −Σ p_i log₂(p_i), p_i = n_i / N, in bits [1][2][3].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Normalised entropy (H / log₂ N) and the entropy-rank ratio of [2]; users should rely on the raw bits value and normalise externally if needed.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| L < k (e.g. `AC`, k=5) | 0 | No k-mers exist; entropy of empty distribution is 0 |
| Single repeated k-mer (`AAAA`, k=2) | 0 | Deterministic distribution, p=1 [3] |
| All-distinct k-mers (`ACGT`, k=2) | log₂(N) | Uniform distribution maximum [3] |
| k < 1 | `ArgumentOutOfRangeException` | Invalid window |
| null DnaSequence | `ArgumentNullException` | Contract |
| null/empty string | 0 | String overload contract |

### 6.2 Limitations

The metric does not normalise by the maximum (log₂ N), so values from sequences of different lengths are not directly comparable; the implementation does not validate the residue alphabet (any character is treated as part of a k-mer). It models only k-mer frequency, not positional structure or reverse-complement equivalence.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double h = SequenceComplexity.CalculateKmerEntropy(new DnaSequence("ATATAT"), k: 2);
// h == 0.9709505944546686  (AT appears 3×, TA 2× among N = 5 dimers)
```

**Numerical walk-through:** `ATATAT`, k=2 → dimers AT,TA,AT,TA,AT ⇒ AT=3, TA=2, N=5. p = 0.6, 0.4. H = −0.6·log₂0.6 − 0.4·log₂0.4 = 0.4421793565 + 0.5287712380 = 0.9709505945 bits.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceComplexity_CalculateKmerEntropy_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/SequenceComplexity_CalculateKmerEntropy_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [SEQ-COMPLEX-KMER-001-Evidence.md](../../../docs/Evidence/SEQ-COMPLEX-KMER-001-Evidence.md)

## 8. References

1. Li, H. 2025. Finding low-complexity DNA sequences with longdust. arXiv:2509.07357. https://arxiv.org/pdf/2509.07357
2. Çakır, et al. 2025. Entropy–Rank Ratio: A Novel Entropy-Based Perspective for DNA Complexity and Classification. arXiv:2511.05300. https://arxiv.org/html/2511.05300
3. Shannon, C. E. 1948. A Mathematical Theory of Communication. Bell System Technical Journal 27. https://en.wikipedia.org/wiki/Entropy_(information_theory)
