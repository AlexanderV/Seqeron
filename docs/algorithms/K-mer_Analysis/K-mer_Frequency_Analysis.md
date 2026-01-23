# K-mer Frequency Analysis

## Overview

K-mer frequency analysis extends basic k-mer counting by computing normalized frequency distributions (where frequencies sum to 1.0), the k-mer spectrum (frequency-of-frequency distribution), and k-mer entropy (Shannon entropy over k-mer distribution). These derived metrics are essential for sequence comparison, genome assembly quality assessment, and metagenomics binning.

## Definition

### K-mer Frequencies (Normalized Counts)

For a sequence with counted k-mers, **k-mer frequencies** are the normalized probabilities:

$$f_i = \frac{c_i}{\sum_j c_j}$$

Where:
- $c_i$ = count of k-mer $i$
- $\sum_j c_j$ = total k-mer count = $L - k + 1$

**Invariant:** $\sum_i f_i = 1.0$

### K-mer Spectrum

The **k-mer spectrum** shows how many k-mers appear with each frequency (multiplicity). It is a histogram mapping count → number of k-mers with that count.

**Example:** For sequence "ACGTACGT" with k=4:
- ACGT appears 2 times
- CGTA, GTAC, TACG each appear 1 time

Spectrum: {1: 3, 2: 1} (3 k-mers appear once, 1 k-mer appears twice)

The k-mer spectrum is widely used in genome assembly to estimate genome size, detect sequencing errors, and identify repeat regions. (Wikipedia: K-mer)

### K-mer Entropy (Shannon Entropy)

**K-mer entropy** measures the diversity of k-mer composition using Shannon's information entropy:

$$H = -\sum_i f_i \log_2(f_i)$$

Where $f_i$ is the frequency of k-mer $i$.

**Properties:**
- **Minimum entropy (H = 0):** All k-mers are identical (homopolymer)
- **Maximum entropy (H = log₂(n)):** All k-mers appear with equal frequency, where n = number of unique k-mers
- For DNA with k=1: Maximum entropy = log₂(4) = 2 bits

**Key principle from Shannon (1948):** "Entropy is maximal when all outcomes are equally likely."

## Edge Cases

| Scenario | Frequencies | Spectrum | Entropy |
|----------|-------------|----------|---------|
| Empty sequence | Empty dictionary | Empty dictionary | 0.0 |
| k > sequence length | Empty dictionary | Empty dictionary | 0.0 |
| Single k-mer possible | {kmer: 1.0} | {1: 1} | 0.0 |
| Homopolymer (AAAA, k=2) | {"AA": 1.0} | {3: 1} | 0.0 |
| All distinct k-mers | Equal frequencies | {1: n} | log₂(n) |

## Invariants

1. **Frequency sum invariant:** Sum of all k-mer frequencies = 1.0
2. **Spectrum total invariant:** Sum of (multiplicity × count) over spectrum = L − k + 1
3. **Entropy bounds:** 0 ≤ H ≤ log₂(unique k-mer count)
4. **Homopolymer entropy:** Homopolymer sequences (e.g., "AAAA") have entropy = 0

## Applications

- **Genome assembly:** K-mer spectrum analysis for genome size estimation (Wikipedia: K-mer)
- **Metagenomics binning:** Tetranucleotide (k=4) frequencies as genomic signatures (Teeling et al., 2004)
- **Sequence comparison:** Alignment-free distance metrics using k-mer frequency profiles
- **Sequencing error detection:** Low-frequency k-mers often represent errors

## Implementation Notes

### Current Implementation

| Method | Class | Description |
|--------|-------|-------------|
| `GetKmerFrequencies(string, k)` | KmerAnalyzer | Returns normalized frequencies (0.0–1.0) |
| `GetKmerSpectrum(string, k)` | KmerAnalyzer | Returns frequency-of-frequency distribution |
| `CalculateKmerEntropy(string, k)` | KmerAnalyzer | Returns Shannon entropy in bits |

### Frequency Calculation

1. Count all k-mers using `CountKmers`
2. Compute total = sum of all counts
3. Divide each count by total
4. Return dictionary of {k-mer: frequency}

### Spectrum Calculation

1. Count all k-mers using `CountKmers`
2. For each count value, count how many k-mers have that count
3. Return dictionary of {count: number_of_kmers}

### Entropy Calculation

Uses base-2 logarithm (log₂) for entropy in bits. Handles the edge case where f=0 by convention: 0 × log(0) = 0.

The implementation rounds to 4 decimal places for numerical stability.

## References

1. Wikipedia. "K-mer." https://en.wikipedia.org/wiki/K-mer
2. Wikipedia. "Entropy (information theory)." https://en.wikipedia.org/wiki/Entropy_(information_theory)
3. Shannon, C.E. (1948). "A Mathematical Theory of Communication." Bell System Technical Journal, 27(3), 379–423.
4. Rosalind. "K-mer Composition." https://rosalind.info/problems/kmer/
5. Teeling, H. et al. (2004). "TETRA: a web-service and a stand-alone program for the analysis and comparison of tetranucleotide usage patterns in DNA sequences." BMC Bioinformatics, 5:163.
6. Chor, B. et al. (2009). "Genomic DNA k-mer spectra: models and modalities." Genome Biology, 10(10): R108.

## Test Unit

- **ID:** KMER-FREQ-001
- **Methods:** GetKmerSpectrum, GetKmerFrequencies, CalculateKmerEntropy
