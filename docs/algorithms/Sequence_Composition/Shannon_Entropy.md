# Shannon Entropy for Biological Sequences

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-ENTROPY-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Shannon entropy measures information content or uncertainty in a sequence. For DNA, it quantifies nucleotide diversity and is commonly used to identify low-complexity regions, repetitive sequence, and conservation patterns. In this repository, the canonical DNA-focused implementation counts only `A/T/G/C`, while a separate general-purpose implementation counts all letters.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Shannon entropy summarizes symbol diversity in a distribution, and in sequence analysis it is used for low-complexity detection, sequence logos, and comparisons between more and less compositionally constrained regions. The original document also highlights k-mer entropy as a higher-order extension that captures local sequence structure beyond single-base composition. Sources: Shannon (1948), Cover & Thomas (1991), Schneider & Stephens (1990), Wikipedia (Entropy, Sequence logo, K-mer).

### 2.2 Core Model

Base-composition entropy is:

$$
H(X) = -\sum_{i=1}^{n} p(x_i) \log_2 p(x_i)
$$

For DNA with 4 bases, the maximum entropy is `log2(4) = 2` bits. The k-mer extension is:

$$
H_k = -\sum_{kmer} p(kmer) \log_2 p(kmer)
$$

where the maximum possible value over a full DNA k-mer alphabet is `log2(4^k) = 2k` bits.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | DNA entropy from the canonical implementation lies in `[0, 2]` bits when only `A/T/G/C` are considered | The implementation counts a 4-symbol alphabet |
| INV-02 | Empty sequences return `0.0` | The source short-circuits before computing frequencies |
| INV-03 | Homopolymers have entropy `0.0` | One symbol has probability `1.0` and all others have probability `0` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | Sequence whose entropy is analyzed | Null `DnaSequence` input throws `ArgumentNullException`; empty string returns `0.0` |
| `k` | `int` | `2` | K-mer size for the `CalculateKmerEntropy(...)` overload | Typed overload throws `ArgumentOutOfRangeException` when `k < 1` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `entropy` | `double` | Shannon entropy in bits |

### 3.3 Preconditions and Validation

`SequenceComplexity.CalculateShannonEntropy(DnaSequence)` throws `ArgumentNullException` for null input. The raw-string overload returns `0.0` for null or empty strings and uppercases the sequence before analysis. `SequenceComplexity.CalculateKmerEntropy(DnaSequence, int)` throws `ArgumentOutOfRangeException` for `k < 1` and returns `0.0` when sequence length is shorter than `k`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the sequence to uppercase.
2. Count base frequencies over the canonical DNA alphabet `A/T/G/C`.
3. Convert counts to probabilities over the counted total.
4. Sum `-p * log2(p)` over non-zero probabilities.
5. For k-mer entropy, count overlapping k-mers and apply the same entropy formula to their frequency distribution.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

DNA-specific entropy values preserved from the original document:

| Distribution | Entropy (bits) | Description |
|-------------|----------------|-------------|
| Uniform `25%` each | `2.0` | Maximum DNA entropy |
| Two bases `50/50` | `1.0` | Intermediate |
| Single base `100%` | `0.0` | Minimum entropy |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateShannonEntropy` | `O(n)` | `O(1)` | Counts the four DNA bases |
| `CalculateKmerEntropy` | `O(n)` | `O(u)` | Builds a k-mer count dictionary before computing entropy |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs), [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceComplexity.CalculateShannonEntropy(...)`: Canonical DNA-oriented entropy implementation.
- `SequenceComplexity.CalculateKmerEntropy(...)`: DNA-oriented k-mer entropy implementation.
- `SequenceStatistics.CalculateShannonEntropy(string)`: General-purpose alternative that counts all letters.

### 5.2 Current Behavior

The canonical `SequenceComplexity` implementation counts only `A/T/G/C` and ignores non-standard bases such as `N` or other ambiguity codes. It uses base-2 logarithms and returns `0.0` if no counted DNA bases are present after filtering. The alternative `SequenceStatistics.CalculateShannonEntropy(...)` counts all letters and is therefore more appropriate for non-DNA alphabets or general text.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Shannon entropy over observed symbol frequencies.
- Base-2 logarithms so entropy is expressed in bits.
- K-mer entropy as a higher-order extension over overlapping substrings.

**Intentionally simplified:**

- The canonical `SequenceComplexity` implementation ignores non-`ATGC` symbols; **consequence:** ambiguity codes do not contribute to the entropy value and the reported range remains DNA-specific.

**Not implemented:**

- General-alphabet entropy in the canonical DNA path; **users should rely on:** `SequenceStatistics.CalculateShannonEntropy(string)` when all letters must be counted.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | `0.0` | No information content |
| Single repeated base | `0.0` | No uncertainty |
| Homopolymer such as `AAAA` | `0.0` | One symbol has probability `1.0` |
| Equal use of all four bases | `2.0` | Maximum uncertainty for DNA |
| Sequence length shorter than `k` for k-mer entropy | `0.0` | No k-mers can be extracted |

### 6.2 Limitations

The canonical implementation is intentionally DNA-focused and ignores non-`ATGC` symbols. It therefore does not provide a full general-alphabet entropy measure in the main API and must be paired with the general-purpose `SequenceStatistics` helper when that behavior is needed.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

- Low-complexity region detection.
- Sequence-logo information content via `2 - H` for DNA positions.
- Comparison of coding and non-coding composition patterns.
- Entropy-like masking strategies such as DUST/SEG-style workflows.

## 8. References

1. Shannon, C.E. (1948). "A Mathematical Theory of Communication". Bell System Technical Journal. 27(3): 379–423.
2. Cover, T.M. & Thomas, J.A. (1991). Elements of Information Theory. Wiley.
3. Schneider, T.D. & Stephens, R.M. (1990). "Sequence Logos: A New Way to Display Consensus Sequences". Nucleic Acids Res. 18(20): 6097–6100.
4. Wikipedia. "Entropy (information theory)". https://en.wikipedia.org/wiki/Entropy_(information_theory)
5. Wikipedia. "Sequence logo". https://en.wikipedia.org/wiki/Sequence_logo
6. Wikipedia. "K-mer". https://en.wikipedia.org/wiki/K-mer
