# Linguistic Complexity (LC)

## Overview

Linguistic Complexity (LC) is a measure of the "vocabulary richness" of a nucleotide sequence, reflecting how many distinct subsequences (subwords) appear in the sequence relative to the maximum possible number.

**Test Unit ID:** SEQ-COMPLEX-001
**Implementation:** `SequenceComplexity.CalculateLinguisticComplexity()`

## Mathematical Definition

### Formula (Summation Variant)

The implementation uses the summation variant as described in Troyanskaya et al. (2002) and Orlov & Potapov (2004):

$$LC = \frac{\sum_{i=1}^{m} V_i}{\sum_{i=1}^{m} V_{max,i}}$$

Where:
- $V_i$ = number of distinct subwords of length $i$ observed in the sequence
- $V_{max,i}$ = maximum possible distinct subwords of length $i$
- $m$ = maximum word length parameter (default: 10)

### Maximum Possible Subwords

For a sequence of length $N$ with alphabet size $K = 4$ (DNA):

$$V_{max,i} = \min(K^i, N - i + 1)$$

This accounts for:
1. **Alphabet constraint:** Maximum $4^i$ possible DNA words of length $i$
2. **Position constraint:** Maximum $N - i + 1$ positions where word can occur

## Properties and Invariants

| Property | Value | Source |
|----------|-------|--------|
| **Range** | $0 \leq LC \leq 1$ | Definition |
| **High complexity** | $LC \approx 1$ | Random-like sequences |
| **Low complexity** | $LC \approx 0$ | Highly repetitive sequences |
| **Empty sequence** | $LC = 0$ | Implementation convention |
| **Time complexity** | $O(n \times k)$ | $n$ = sequence length, $k$ = max word length |

## Biological Significance

Low-complexity regions in genomes correlate with:
- Simple sequence repeats (SSRs) / microsatellites
- Tandem repeats
- Imperfect direct and inverted repeats
- Palindrome-hairpin structures
- Potential regulatory sites

High-complexity regions correlate with:
- Coding sequences (exons)
- Functional protein-coding genes

## Implementation Notes

### Current Implementation

The `SequenceComplexity.CalculateLinguisticComplexity()` method:
1. Iterates word lengths from 1 to `min(maxWordLength, sequence.Length)`
2. Uses `HashSet<string>` to count distinct subwords at each length
3. Sums observed vs. maximum possible across all word lengths
4. Returns ratio (0 to 1)

### Edge Cases

| Input | Expected Output | Rationale |
|-------|-----------------|-----------|
| Empty string | 0 | No vocabulary possible |
| Null DnaSequence | `ArgumentNullException` | Guard clause |
| `maxWordLength < 1` | `ArgumentOutOfRangeException` | Invalid parameter |
| Single nucleotide "A" | > 0 | Has vocabulary of size 1 |
| Homopolymer "AAAA..." | Low value | Only 1 distinct N-mer per length |
| Random sequence | High value (~0.8-1.0) | Maximum vocabulary diversity |

## Related Methods in Test Unit

| Method | Type | Description |
|--------|------|-------------|
| `CalculateLinguisticComplexity(DnaSequence, int)` | Canonical | Primary implementation |
| `CalculateLinguisticComplexity(string, int)` | String overload | Convenience wrapper |
| `FindLowComplexityRegions(...)` | Region detection | Uses LC internally |
| `MaskLowComplexity(...)` | Masking | Uses DUST score (related algorithm) |

## References

1. **Trifonov, E.N.** (1990). "Making sense of the human genome." Structure and Methods, Vol. 1. Human Genome Initiative and DNA Recombination. Adenine Press, pp. 69–77.
   - *Original introduction of linguistic complexity concept*

2. **Troyanskaya, O.G., Arbell, O., Koren, Y., Landau, G.M., Bolshoy, A.** (2002). "Sequence complexity profiles of prokaryotic genomic sequences: A fast algorithm for calculating linguistic complexity." Bioinformatics, 18(5), 679–688. DOI: 10.1093/bioinformatics/18.5.679
   - *Fast suffix tree algorithm, summation formula variant*

3. **Orlov, Y.L., Potapov, V.N.** (2004). "Complexity: an internet resource for analysis of DNA sequence complexity." Nucleic Acids Research, 32(Web Server issue), W628–W633. DOI: 10.1093/nar/gkh466
   - *Comparative analysis of complexity measures, implementation details*

4. **Gabrielian, A., Bolshoy, A.** (1999). "Sequence complexity and DNA curvature." Computers & Chemistry, 23(3–4), 263–274.
   - *Extended vocabulary usage measures*

5. **Wikipedia** - "Linguistic sequence complexity"
   - *Summary of approaches and formulas*
