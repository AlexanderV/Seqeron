# Shannon Entropy for Biological Sequences

## Overview

Shannon entropy is a measure of information content (uncertainty) in a sequence. For DNA sequences, it quantifies nucleotide diversity and is commonly used to detect low-complexity regions, repetitive sequences, and conservation patterns.

## Mathematical Definition

Shannon entropy H is defined as:

$$H(X) = -\sum_{i=1}^{n} p(x_i) \log_2 p(x_i)$$

Where:
- $p(x_i)$ is the probability (frequency) of symbol $x_i$
- $n$ is the number of distinct symbols (4 for DNA: A, T, G, C)
- The logarithm base 2 gives entropy in **bits**

**Source:** Shannon, C.E. (1948). "A Mathematical Theory of Communication". Bell System Technical Journal. 27(3): 379–423.

## DNA-Specific Properties

For DNA sequences with 4 nucleotides (A, T, G, C):

| Distribution | Entropy (bits) | Description |
|-------------|----------------|-------------|
| Uniform (25% each) | 2.0 | Maximum entropy |
| Two bases (50% each) | 1.0 | Intermediate |
| Single base (100%) | 0.0 | Minimum entropy (homopolymer) |

**Key invariant:** For DNA, $0 \leq H \leq 2$ bits per base.

**Source:** Wikipedia, "Entropy (information theory)" - Maximum entropy = $\log_2(n)$ where n = alphabet size.

## K-mer Entropy

An extension of Shannon entropy using k-mers (substrings of length k) instead of single bases:

$$H_k = -\sum_{kmer} p(kmer) \log_2 p(kmer)$$

Where:
- K-mers are overlapping substrings of length k
- Maximum possible entropy = $\log_2(4^k) = 2k$ bits

K-mer entropy captures higher-order sequence structure and correlations between adjacent bases that simple base-composition entropy misses.

**Source:** Wikipedia, "K-mer" - k-mers are substrings of length k used in bioinformatics for sequence analysis.

## Applications in Bioinformatics

1. **Low-complexity region detection**: Regions with entropy < threshold indicate repeats/simple sequences
2. **Sequence logos**: Information content at each position = $R_i = \log_2(4) - H_i = 2 - H_i$
3. **Sequence comparison**: Entropy profiles can distinguish coding vs non-coding regions
4. **DUST/SEG algorithms**: Use entropy-like measures for masking repetitive regions

**Source:** Schneider & Stephens (1990). "Sequence Logos: A New Way to Display Consensus Sequences". Nucleic Acids Res. 18(20): 6097–6100.

## Implementation Notes

### Current Implementation: `SequenceComplexity.CalculateShannonEntropy`

The canonical implementation:
- Counts only standard DNA bases (A, T, G, C)
- Ignores non-standard bases (N, ambiguous codes)
- Uses log base 2 for bits
- Returns 0 for empty sequences
- Case-insensitive (converts to uppercase internally)

### Current Implementation: `SequenceComplexity.CalculateKmerEntropy`

The k-mer entropy variant:
- Counts overlapping k-mers
- Returns 0 if sequence length < k
- Requires k ≥ 1

### Alternative Implementation: `SequenceStatistics.CalculateShannonEntropy`

A more general implementation:
- Counts all letters (not just ATGC)
- Suitable for protein sequences or general text
- Different invariant bounds depending on alphabet

**ASSUMPTION:** The `SequenceComplexity` implementation is canonical for DNA analysis; `SequenceStatistics` version is a general-purpose utility.

## Edge Cases

| Case | Expected Result | Rationale |
|------|-----------------|-----------|
| Empty sequence | 0.0 | No information content |
| Single base | 0.0 | No uncertainty |
| Homopolymer (e.g., "AAAA") | 0.0 | Zero diversity |
| All 4 bases equal | 2.0 | Maximum uncertainty |
| Two bases 50/50 | 1.0 | Binary entropy |
| Three bases equal, one absent | ~1.58 | $\log_2(3)$ |
| Sequence length < k (for k-mer) | 0.0 | No k-mers extractable |

## References

1. Shannon, C.E. (1948). "A Mathematical Theory of Communication". Bell System Technical Journal. 27(3): 379–423.
2. Cover, T.M. & Thomas, J.A. (1991). Elements of Information Theory. Wiley.
3. Schneider, T.D. & Stephens, R.M. (1990). "Sequence Logos: A New Way to Display Consensus Sequences". Nucleic Acids Res. 18(20): 6097–6100.
4. Wikipedia. "Entropy (information theory)". https://en.wikipedia.org/wiki/Entropy_(information_theory)
5. Wikipedia. "Sequence logo". https://en.wikipedia.org/wiki/Sequence_logo
6. Wikipedia. "K-mer". https://en.wikipedia.org/wiki/K-mer
