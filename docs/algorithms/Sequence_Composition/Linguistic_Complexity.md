# Linguistic Complexity (LC)

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-COMPLEX-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Linguistic complexity measures how many distinct subsequences appear in a sequence relative to how many could appear in principle. In this repository, the implementation computes the summation variant over word lengths from 1 up to a configurable maximum and uses it as a DNA-oriented complexity metric for low-complexity analysis. The current implementation follows the definition directly with hash-based subword enumeration rather than the suffix-tree optimization discussed in some of the cited literature.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Linguistic complexity is a vocabulary-richness measure for nucleotide sequences. Low values are associated with simple repeats, tandem repeats, and other low-complexity patterns, while higher values are associated with more diverse sequence content such as coding regions. Sources: Trifonov (1990), Troyanskaya et al. (2002), Orlov & Potapov (2004), Gabrielian & Bolshoy (1999), Wikipedia (Linguistic sequence complexity).

### 2.2 Core Model

The implementation uses the summation variant:

$$
LC = \frac{\sum_{i=1}^{m} V_i}{\sum_{i=1}^{m} V_{max,i}}
$$

where `V_i` is the number of distinct observed subwords of length `i`, `V_{max,i}` is the maximum possible number of distinct subwords of that length, and `m` is the maximum word length parameter. For DNA with alphabet size `K = 4`:

$$
V_{max,i} = \min(K^i, N - i + 1)
$$

where `N` is sequence length.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 <= LC <= 1` for DNA-alphabet inputs | Observed distinct counts cannot exceed the DNA-theoretical maximum when the input alphabet matches the hard-coded `K = 4` denominator |
| INV-02 | Empty sequences return `0` | The implementation short-circuits before accumulating counts |
| INV-03 | Word lengths are limited to `min(maxWordLength, sequence.Length)` | The source explicitly caps the loop bound |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA-oriented sequence to analyze | Null `DnaSequence` input throws `ArgumentNullException`; empty string returns `0` |
| `maxWordLength` | `int` | `10` | Maximum subword length included in the summation | `DnaSequence` overload throws `ArgumentOutOfRangeException` when `< 1` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `lc` | `double` | Linguistic-complexity ratio; for DNA-alphabet inputs it lies between `0` and `1` |

### 3.3 Preconditions and Validation

`CalculateLinguisticComplexity(DnaSequence, int)` throws `ArgumentNullException` for null sequences and `ArgumentOutOfRangeException` when `maxWordLength < 1`. The raw-string overload returns `0` for null or empty input and uppercases the sequence before analysis, but it does not validate or filter the alphabet beyond that normalization.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. For each word length from 1 to `min(maxWordLength, sequence.Length)`, enumerate all overlapping subwords.
3. Count distinct subwords for that length with a `HashSet<string>`.
4. Compute the maximum possible count for that length using `min(4^i, N - i + 1)`.
5. Sum observed and possible counts across all lengths and return their ratio.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateLinguisticComplexity` | `O(n × k^2)` effective | `O(u)` | `n` is sequence length, `k` is the effective maximum word length, and the direct implementation allocates and hashes substrings of lengths `1..k` for each window |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateLinguisticComplexity(DnaSequence, int)`: Canonical typed overload.
- `SequenceComplexity.CalculateLinguisticComplexity(string, int)`: Raw-string overload.
- `SequenceComplexity.FindLowComplexityRegions(...)`: Uses complexity metrics downstream.
- `SequenceComplexity.MaskLowComplexity(...)`: Related masking workflow using DUST score.

### 5.2 Current Behavior

The current implementation counts distinct subwords with `HashSet<string>` collections for each word length and sums the observed and possible totals directly. It allocates and hashes fresh substrings for each tested window length rather than using a suffix-tree index. The typed overload enforces `maxWordLength >= 1`, while the raw-string overload uppercases input and delegates to the same core computation without alphabet validation. The denominator remains hard-coded to the DNA alphabet size `4`, so raw-string inputs containing other symbols can exceed the DNA-bounded `[0, 1]` interpretation. The effective word-length range is capped at sequence length.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Summation-form linguistic complexity over subword lengths.
- Maximum distinct-subword bounds based on both alphabet size and positional availability.
- DNA-oriented complexity scoring in the range `[0, 1]`.

**Intentionally simplified:**

- The implementation assumes a DNA alphabet of size 4; **consequence:** the metric is not generalized to arbitrary alphabets in this code path.
- The raw-string overload accepts arbitrary uppercase symbols while still using the DNA denominator `4^i`; **consequence:** callers should treat the reported value as DNA-oriented and not assume the usual `[0, 1]` bound for non-ACGT inputs.

**Not implemented:**

- The fast suffix-tree algorithm described in some cited literature; **users should rely on:** the current direct enumeration path, which prioritizes clarity over that optimization.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Current source uses hash-based subword enumeration rather than the fast suffix-tree approach highlighted by Troyanskaya et al. (2002) | Deviation | Runtime follows direct enumeration rather than a suffix-tree optimization | accepted | Confirmed from `SequenceComplexity.CalculateLinguisticComplexityCore(...)` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty string | Returns `0` | No subword vocabulary exists |
| Null `DnaSequence` | Throws `ArgumentNullException` | Explicit guard |
| `maxWordLength < 1` on typed overload | Throws `ArgumentOutOfRangeException` | Explicit validation |
| Single nucleotide such as `A` | Returns a positive value | One distinct 1-mer exists |
| Homopolymer sequence | Returns a low value | Only one word per length is observed |
| Random-like sequence | Returns a high value | Observed vocabulary approaches the maximum |
| Raw-string input with non-ACGT symbols | May exceed the usual DNA-bounded interpretation | Observed words can include symbols outside the hard-coded DNA denominator |

### 6.2 Limitations

The current implementation is DNA-specific and uses direct `HashSet<string>` enumeration rather than the faster suffix-tree approach discussed in some of the cited papers. Runtime therefore includes repeated substring allocation and hashing across the tested word lengths, and memory usage grows with the number of distinct observed subwords. The raw-string overload also accepts arbitrary uppercase symbols without reconciling the denominator to a larger alphabet.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

- Low-complexity region detection in repetitive DNA.
- Characterization of microsatellites, tandem repeats, and palindrome-hairpin-rich segments.
- Contrast between low-complexity and coding-like sequence regions.

## 8. References

1. Trifonov, E.N. (1990). "Making sense of the human genome."
2. Troyanskaya, O.G., Arbell, O., Koren, Y., Landau, G.M., Bolshoy, A. (2002). "Sequence complexity profiles of prokaryotic genomic sequences: A fast algorithm for calculating linguistic complexity." Bioinformatics, 18(5), 679–688.
3. Orlov, Y.L., Potapov, V.N. (2004). "Complexity: an internet resource for analysis of DNA sequence complexity." Nucleic Acids Research, 32(Web Server issue), W628–W633.
4. Gabrielian, A., Bolshoy, A. (1999). "Sequence complexity and DNA curvature." Computers & Chemistry, 23(3–4), 263–274.
5. Wikipedia - "Linguistic sequence complexity".
   - *Summary of approaches and formulas*
