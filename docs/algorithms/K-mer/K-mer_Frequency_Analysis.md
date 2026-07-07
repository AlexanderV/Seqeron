# K-mer Frequency Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer Analysis |
| Test Unit ID | KMER-FREQ-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

K-mer frequency analysis extends basic k-mer counting by deriving normalized k-mer frequencies, the k-mer spectrum, and k-mer entropy. These quantities are useful for sequence comparison, genome-assembly quality assessment, and metagenomics signatures. In this repository, all three metrics are built directly from exact k-mer counts returned by `KmerAnalyzer.CountKmers(...)`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Normalized k-mer frequencies convert counts into probabilities, the k-mer spectrum records how many k-mers occur with each multiplicity, and Shannon entropy summarizes the diversity of the resulting distribution. The original document also notes the use of k-mer spectra for assembly and error detection and tetranucleotide frequencies for metagenomics signatures. Sources: Wikipedia (K-mer, Entropy), Shannon (1948), Teeling et al. (2004), Chor et al. (2009), Rosalind.

### 2.2 Core Model

Normalized frequency for k-mer `i` is:

$$
f_i = \frac{c_i}{\sum_j c_j}
$$

where `c_i` is the observed count of k-mer `i`. The k-mer spectrum is the histogram mapping `count -> number of k-mers with that count`. Shannon k-mer entropy is:

$$
H = -\sum_i f_i \log_2(f_i)
$$

with the convention that terms with `f_i = 0` contribute `0`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The sum of all returned frequencies is `1.0` when at least one k-mer exists | Frequencies are each divided by the total count |
| INV-02 | The spectrum total satisfies `Σ(count × multiplicity) = L - k + 1` when `k <= L` | Spectrum bins are derived from exact k-mer counts |
| INV-03 | `0 <= H <= log2(unique k-mer count)` | Entropy is computed from a discrete probability distribution |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Sequence whose k-mer distribution is analyzed | Null or empty string yields empty outputs or zero entropy |
| `k` | `int` | required | K-mer length | `k <= 0` throws through the underlying count routine |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `frequencies` | `Dictionary<string, double>` | Normalized frequency per observed k-mer |
| `spectrum` | `Dictionary<int, int>` | Histogram mapping count to number of k-mers |
| `entropy` | `double` | Shannon entropy in bits |

### 3.3 Preconditions and Validation

All three metrics delegate to `CountKmers(...)` for input handling. Null or empty sequences yield empty dictionaries and entropy `0.0`. If `k` exceeds sequence length, the count dictionary is empty and entropy is `0.0`. If `k <= 0`, the underlying counting routine throws `ArgumentOutOfRangeException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count all k-mers in the sequence.
2. Compute the total count and divide each count by that total to obtain normalized frequencies.
3. Invert the count dictionary to build the multiplicity spectrum.
4. Sum `-f * log2(f)` over the non-zero frequencies to obtain entropy.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GetKmerFrequencies` | `O(n)` | `O(u)` | Derived from exact counts |
| `GetKmerSpectrum` | `O(n)` | `O(u)` | Iterates over the count values |
| `CalculateKmerEntropy` | `O(n)` | `O(u)` | Builds on normalized frequencies |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.GetKmerFrequencies(string, int)`: Returns normalized frequencies in `[0.0, 1.0]`.
- `KmerAnalyzer.GetKmerSpectrum(string, int)`: Returns the count-of-counts histogram.
- `KmerAnalyzer.CalculateKmerEntropy(string, int)`: Returns Shannon entropy in bits.

### 5.2 Current Behavior

The current implementation always computes these metrics from exact k-mer counts. Frequency normalization uses the sum of observed counts, not the theoretical number of possible k-mers. Entropy uses `Math.Log2` and skips zero-frequency terms by iterating only over observed frequencies.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Frequency normalization by total observed k-mer count.
- Spectrum construction as a histogram of k-mer multiplicities.
- Shannon entropy over the observed k-mer distribution using base-2 logarithms.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The original document described 4-decimal entropy rounding for numerical stability, but the current source returns the raw double sum without an explicit rounding step | Deviation | Reported entropy may include full floating-point precision | accepted | Confirmed from `CalculateKmerEntropy(...)` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Empty frequency map, empty spectrum, entropy `0.0` | No k-mers exist |
| `k > sequence.Length` | Empty frequency map, empty spectrum, entropy `0.0` | No valid windows exist |
| Single possible k-mer | Frequency `1.0`, spectrum `{1: 1}`, entropy `0.0` | The distribution has one outcome |
| Homopolymer such as `AAAA`, `k = 2` | Frequency `{"AA": 1.0}`, spectrum `{3: 1}`, entropy `0.0` | All windows are identical |

### 6.2 Limitations

The current implementation analyzes only observed k-mers and does not smooth the distribution or normalize against theoretical k-mer space. As with the underlying counting routine, memory usage grows with the number of unique observed k-mers.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

- Genome assembly through k-mer spectrum analysis.
- Metagenomics binning via tetranucleotide signatures.
- Alignment-free sequence comparison using k-mer profiles.
- Sequencing-error detection from low-frequency k-mers.

## 8. References

1. Wikipedia. "K-mer." https://en.wikipedia.org/wiki/K-mer
2. Wikipedia. "Entropy (information theory)." https://en.wikipedia.org/wiki/Entropy_(information_theory)
3. Shannon, C.E. (1948). "A Mathematical Theory of Communication." Bell System Technical Journal, 27(3), 379–423.
4. Rosalind. "K-mer Composition." https://rosalind.info/problems/kmer/
5. Teeling, H. et al. (2004). "TETRA: a web-service and a stand-alone program for the analysis and comparison of tetranucleotide usage patterns in DNA sequences." BMC Bioinformatics, 5:163.
6. Chor, B. et al. (2009). "Genomic DNA k-mer spectra: models and modalities." Genome Biology, 10(10): R108.
