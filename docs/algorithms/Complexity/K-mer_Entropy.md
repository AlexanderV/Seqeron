# K-mer Entropy

| Field | Value |
|-------|-------|
| Algorithm Group | Complexity |
| Test Unit ID | SEQ-COMPLEX-KMER-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

K-mer entropy measures the sequence complexity of a DNA string as the Shannon entropy (in bits) of the frequency distribution of its overlapping k-mers. It quantifies how uniformly the k-mers are distributed: a low value indicates a few dominant k-mers (repeats, homopolymers, low complexity), while a high value indicates a near-uniform k-mer distribution (high complexity) [1]. The computation is exact and deterministic for a given sequence and k.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Detecting low-complexity DNA regions is a standard pre-processing step for sequence alignment and search. The Shannon entropy of overlapping n-word (k-mer) frequencies — the *block entropy* H_n — is the standard information-theoretic complexity measure of DNA [4][5]; read filters such as BBDuk compute exactly this quantity (then divide by ln N) [6]. Note that `longdust` [1] uses a different, Poisson composite-likelihood score (Σ log c(t)! − f(ℓ/4^k)), not Shannon entropy; it is cited here only for the overlapping-k-mer count ℓ = L − k + 1. The entropy saturates at its maximum for random uniform sequences and drops toward zero for repetitive ones [2].

### 2.2 Core Model

Decompose a sequence of length L into its overlapping k-mers using a sliding window of step 1, giving N = L − k + 1 k-mers [1][6]. Let n_i be the count of the i-th distinct k-mer and p_i = n_i / N its relative frequency, so Σ p_i = 1. The Shannon entropy is

> H = − Σ_i p_i · log₂(p_i)   (bits) [3][4]

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
- `SequenceComplexity.CalculateKmerEntropyCore(string, int)` (private): counts overlapping k-mers with `KmerAnalyzer.CountKmers` and applies the shared entropy kernel `ShannonEntropyBits`, which delegates to the canonical `StatisticsHelper.ShannonIndex` (natural log) ÷ ln 2 — scipy `entropy(counts, base=2)`'s computation (2000 random cases vs scipy: max |Δ| 3.6e-14, summation-order rounding).

### 5.2 Current Behavior

K-mers are enumerated by the canonical `KmerAnalyzer.CountKmers` (single linear scan, `Dictionary<string,int>`); entropy is then computed over the count values by the private `ShannonEntropyBits` kernel. The repository suffix tree was evaluated and **not** used: this is a single linear pass building a full frequency table (every position visited once), not a repeated occurrence-query workload, so a suffix tree adds construction overhead without changing the linear cost or the required output.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Overlapping k-mer decomposition with N = L − k + 1 [1].
- H = −Σ p_i log₂(p_i), p_i = n_i / N, in bits [3][4][6] — the plug-in (maximum-likelihood) estimate.
- K-mer counts come from the canonical counter `KmerAnalyzer.CountKmers` (KMER-COUNT-001); the entropy kernel is shared with `CalculateShannonEntropy`.

**Intentionally simplified:**

- (none)

**Not implemented:**

- Finite-sample bias correction of block entropies [5] (the method returns the raw plug-in value, which underestimates the source entropy when N ≪ 4^k).
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
4. Herzel, H., Ebeling, W., Schmitt, A. O. 1994. Entropies of biosequences: the role of repeats. Phys. Rev. E 50:5061–5071. https://doi.org/10.1103/PhysRevE.50.5061
5. Schmitt, A. O., Herzel, H. 1997. Estimating the entropy of DNA sequences. J. Theor. Biol. 188:369–377. https://doi.org/10.1006/jtbi.1997.0493
6. Bushnell, B. BBMap/BBDuk `EntropyTracker.java` (reference implementation; pk = count/(window − k + 1)). https://github.com/BioInfoTools/BBMap/blob/master/current/structures/EntropyTracker.java
