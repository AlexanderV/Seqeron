# Linguistic Complexity (LC)

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-COMPLEX-001 |
| Related Projects | N/A |
| Implementation Status | Complete |
| Last Reviewed | 2026-09-30 |

## 1. Overview

Linguistic complexity measures how many distinct subsequences appear in a sequence relative to how many could appear in principle. In this repository, the implementation computes the summation variant over word lengths from 1 up to a configurable maximum and uses it as a DNA-oriented complexity metric for low-complexity analysis. Small word-length limits use direct hash-based subword enumeration; larger limits (m > 12, including Troyanskaya's all-length LC) count distinct subwords from the suffix tree in linear time, as in Troyanskaya et al. (2002). Both paths give identical values.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Linguistic complexity is a vocabulary-richness measure for nucleotide sequences. Low values are associated with simple repeats, tandem repeats, and other low-complexity patterns, while higher values are associated with more diverse sequence content such as coding regions. Sources: Trifonov (1990), Troyanskaya et al. (2002), Orlov & Potapov (2004), Gabrielian & Bolshoy (1999), Wikipedia (Linguistic sequence complexity).

### 2.2 Core Model

The implementation uses the summation variant:

$$
LC = \frac{\sum_{i=1}^{m} V_i}{\sum_{i=1}^{m} V_{max,i}}
$$

where `V_i` is the number of distinct observed subwords of length `i`, `V_{max,i}` is the maximum possible number of distinct subwords of that length, and `m` is the maximum word length parameter. With alphabet size `K` (Troyanskaya et al. 2002; Rosalind LING: "for an alphabet of size a"):

$$
V_{max,i} = \min(K^i, N - i + 1)
$$

where `N` is sequence length. The implementation takes `K = |{A, C, G, T} ∪ symbols(s)|` (U replaces T when U occurs and T does not): pure DNA/RNA gives `K = 4`; every other symbol present (N, IUPAC codes, gaps, …) enlarges the alphabet, so `V_i ≤ V_max,i` and `LC ≤ 1` for every input (e.g. `ACGTN` → 15/15 = 1.0, `ATGCATGCNN` → 44/50 = 0.88).

This word-length-limited sum is the Orlov & Potapov (2004) CL (`m ≤ N`); with `m ≥ N` it is exactly the Troyanskaya et al. (2002) LC `A(s)/M(s)` over all lengths (Rosalind LING sample `ATTTGGATT` → 0.875). It is distinct from Trifonov's (1990) product form `C = Π U_i` (implemented e.g. by the R package universalmotif, method "Trifonov"); the two are not interchangeable.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 <= LC <= 1` for every input | Every observed symbol is in the alphabet of size `K`, so `V_i <= min(K^i, N-i+1)` |
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
| `lc` | `double` | Linguistic-complexity ratio in `[0, 1]` |

### 3.3 Preconditions and Validation

`CalculateLinguisticComplexity(DnaSequence, int)` throws `ArgumentNullException` for null sequences and `ArgumentOutOfRangeException` when `maxWordLength < 1`. The raw-string overload returns `0` for null or empty input and uppercases the sequence before analysis, but it does not validate or filter the alphabet beyond that normalization.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. Let `m = min(maxWordLength, sequence.Length)`. If `m ≤ 12`, count the distinct overlapping subwords of each length as the key count of the canonical `KmerAnalyzer.CountKmers` tally.
3. Otherwise build (or reuse the cached `DnaSequence.SuffixTree`) suffix tree and call the shared `ISuffixTree.CountDistinctSubstringsByLength(m)` (SuffixTree project): for every edge spanning depths `d+1..d+len` (leaf edges excluding the terminator), add 1 to `V_i` for each covered `i ≤ m` (difference array).
4. Compute the maximum possible count for that length using `min(K^i, N - i + 1)` (K = 4 for DNA/RNA, extended by any other symbol present).
5. Sum observed and possible counts across all lengths and return their ratio.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateLinguisticComplexity` (m ≤ 12) | `O(n × m^2)` | `O(n × m)` | direct substring hashing |
| `CalculateLinguisticComplexity` (m > 12) | `O(n)` (+ tree build) | `O(n)` | suffix-tree edge-depth counting (Troyanskaya 2002) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceComplexity.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs)

- `SequenceComplexity.CalculateLinguisticComplexity(DnaSequence, int)`: Canonical typed overload.
- `SequenceComplexity.CalculateLinguisticComplexity(string, int)`: Raw-string overload.
- `SequenceComplexity.FindLowComplexityRegions(...)`: Uses complexity metrics downstream.
- `SequenceComplexity.MaskLowComplexity(...)`: Related masking workflow using DUST score.

### 5.2 Current Behavior

For `m ≤ 12` the implementation counts distinct subwords as the key count of the canonical `KmerAnalyzer.CountKmers` tally; for larger `m` it counts them from the suffix tree in linear time (terminator excluded). `V_max,i` uses a saturating `K^i` (multiplication stops once it exceeds `N`, so no overflow for any `K ≤ 65536` or `i`). The typed overload enforces `maxWordLength >= 1` (and only admits ACGT, so `K = 4`), while the raw-string overload uppercases input and uses `K = |{A,C,G,T/U} ∪ symbols|` (2026-09-30, B04 F20: previously `K = 4` always, so `ACGTN` gave 15/14 > 1). The effective word-length range is capped at sequence length.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Summation-form linguistic complexity over subword lengths.
- Maximum distinct-subword bounds based on both alphabet size and positional availability.
- Complexity scoring in the range `[0, 1]` for any input alphabet (DNA/RNA `K = 4`).

**Intentionally simplified:**

- None for the alphabet: non-ACGT symbols extend the alphabet size `K` (not filtered out).

- The linear-time suffix-tree counting of Troyanskaya et al. (2002) for `m > 12` (incl. all-length LC).

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Hash enumeration is used for `m ≤ 12`, suffix tree above | Implementation choice | None on values (both exact; locked by LC-14..LC-17) | resolved 2026-09 | Previously the suffix-tree path was missing (O(N²) for all-length LC) |

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
| Raw-string input with non-ACGT symbols | Alphabet enlarged (`ACGTN` → 1.0, `ATGCATGCNN` → 0.88, `NNNNNNNN` → 8/33); never > 1 | `K = |{A,C,G,T/U} ∪ symbols|` |

### 6.2 Limitations

The alphabet is inferred from the sequence (nucleotide alphabet ∪ observed symbols); a caller wanting a different fixed alphabet (e.g. 20 amino acids for a short peptide lacking some residues) cannot pass it explicitly.

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
6. Rosalind, problem LING "Linguistic Complexity of a Genome" (sample ATTTGGATT → 0.875).
7. universalmotif (R), `R/sequence_complexity.R` — Trifonov product-form reference implementation.
   - *Summary of approaches and formulas*
