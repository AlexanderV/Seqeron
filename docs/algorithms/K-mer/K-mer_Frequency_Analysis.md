# K-mer Frequency Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer Analysis |
| Test Unit ID | KMER-FREQ-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-09-28 |

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

where `c_i` is the observed count of k-mer `i`; since every one of the `L - k + 1` overlapping windows is counted, `Σ c_j = L - k + 1` (same denominator as scikit-bio `kmer_frequencies(k, overlap=True, relative=True)`). The k-mer spectrum is the histogram mapping `count -> number of distinct k-mers with that count` [6] (Jellyfish `histo` semantics; only non-zero bins are returned, no upper cap bin, single strand). The full `jellyfish histo` contract — `-l/--low` (default 1), `-h/--high` (default 10000), `-i/--increment` (default 1), the catch-all cap bin and `-f/--full` — is available as `GetKmerHistogram` (§5.1, §7.3) [7]. Shannon k-mer entropy is:

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
| `k` | `int` | required | K-mer length | `k <= 0` throws through the underlying count routine for non-empty input (null/empty input short-circuits to empty/0) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `frequencies` | `Dictionary<string, double>` | Normalized frequency per observed k-mer |
| `spectrum` | `Dictionary<int, int>` | Histogram mapping count to number of k-mers |
| `entropy` | `double` | Shannon entropy in bits |

### 3.3 Preconditions and Validation

All three metrics delegate to `CountKmers(...)` for input handling. Null or empty sequences yield empty dictionaries and entropy `0.0`. If `k` exceeds sequence length, the count dictionary is empty and entropy is `0.0`. If `k <= 0` and the sequence is non-empty, the underlying counting routine throws `ArgumentOutOfRangeException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count all k-mers in the sequence.
2. Compute the total count and divide each count by that total to obtain normalized frequencies.
3. Invert the count dictionary to build the multiplicity spectrum.
4. Sum `-f * log2(f)` over the non-zero frequencies to obtain entropy.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GetKmerFrequencies` | `O(n·k)` | `O(u·k)` | Derived from exact counts; each window builds/hashes a k-length string |
| `GetKmerSpectrum` | `O(n·k)` | `O(u·k)` | Iterates over the count values |
| `CalculateKmerEntropy` | `O(n·k)` | `O(u·k)` | Delegates to the canonical `SequenceComplexity.CalculateKmerEntropy` (SEQ-COMPLEX-KMER-001; `StatisticsHelper.ShannonIndex` ÷ ln 2), bit-identical; own contract kept for empty input / k ≤ 0 (B06 KMER-STATS-001) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.GetKmerFrequencies(string, int)`: Returns normalized frequencies in `[0.0, 1.0]`.
- `KmerAnalyzer.GetKmerSpectrum(string, int)`: Returns the count-of-counts histogram.
- `KmerAnalyzer.GetKmerSpectrum(string, int, KmerCountingOptions)`: The same over literal / ACGT-only / canonical (`count -C`) counts.
- `KmerAnalyzer.GetKmerHistogram(string, int, KmerCountingOptions = default, long low = 1, long high = 10000, long increment = 1, bool full = false)` and `GetKmerHistogram(IEnumerable<int> kmerCounts, …)`: `jellyfish count [-C]` + `jellyfish histo -l -h -i [-f]`, returning ordered `KmerHistogramBin(Bin, Frequency)` rows — exactly the lines Jellyfish prints (B06 audit round 1 WP3, F12).
- `KmerAnalyzer.CalculateKmerEntropy(string, int)`: Returns Shannon entropy in bits.

### 5.2 Current Behavior

The current implementation always computes these metrics from exact k-mer counts. Frequency normalization uses the sum of observed counts, not the theoretical number of possible k-mers. Entropy delegates to `SequenceComplexity.CalculateKmerEntropy` (canonical `StatisticsHelper.ShannonIndex` in nats ÷ ln 2, the `scipy.stats.entropy(counts, base=2)` computation) over the observed counts; null/empty input returns 0 for any k and k ≤ 0 throws only for non-empty input. Until B06 (KMER-STATS-001) it was a separate `Math.Log2` loop over `GetKmerFrequencies`; the two differ by ≤ 4.6e-14 bits (20 000 random tables).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Frequency normalization by total observed k-mer count.
- Spectrum construction as a histogram of k-mer multiplicities.
- Jellyfish `histo` binning (`sub_commands/histo_main.cc`, verbatim): `base = inc >= low ? 0 : low − inc`, `ceil = high + inc`, `nb_buckets = (ceil + inc − base) / inc`; a count `< base` goes to bucket 0, `> ceil` to the last bucket, else to `(count − base) / inc`; bucket i is labelled `base + i·inc`; zero rows only with `--full`. Hence the last bucket (label ≥ high) is the cap for every count above `high`, counts below `low` are pooled in the first bucket, and with `inc ≥ low` the first label is 0. `high < low` is rejected (Jellyfish: "High count value must be >= to low count value"); `inc = 0` (a division by zero in Jellyfish) is rejected.
- Shannon entropy over the observed k-mer distribution using base-2 logarithms.

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The original document described 4-decimal entropy rounding for numerical stability, but the current source returns the raw double sum without an explicit rounding step | Deviation | Reported entropy may include full floating-point precision | accepted | Confirmed from `CalculateKmerEntropy(...)`; matches `scipy.stats.entropy(counts, base=2)` to 1e-12 (review 2026-09) |

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

### 7.3 Jellyfish `histo` cross-check (B06 audit round 1, WP3)

Reference: the real **Jellyfish 2.3.1** binary (`apt jellyfish 2.3.1-3build1`): `jellyfish count -m k -s 10000 -t 1 [-C]` then `jellyfish histo [-l] [-h] [-i] [-f]`. 3 inputs (Rosalind KMER sample k=4, BA1B sample k=4, `A^31 CGTACGTACGTACGTTTGCA` k=3) × {plain, `-C`} × 12 option sets = 72 runs; a Python replica of `histo_main.cc` reproduces all 72, and `GetKmerHistogram` equals every row (`KmerAnalyzer_HistogramClumpWindowsFilters_Tests`). Selected rows (bin frequency …):

| Input | Mode | histo options | Jellyfish 2.3.1 output |
|---|---|---|---|
| Rosalind k=4 | `-C` | (defaults) | 1 23 2 34 3 27 4 25 5 7 6 5 7 4 9 3 10 2 |
| Rosalind k=4 | `-C` | `-h 5` | 1 23 2 34 3 27 4 25 5 7 6 14 (cap bin 6 = counts ≥ 6) |
| Rosalind k=4 | `-C` | `-i 2` | 0 23 2 61 4 32 6 9 8 3 10 2 |
| Rosalind k=4 | `-C` | `-l 3 -h 8 -i 2` | 1 57 3 52 5 12 7 4 9 5 |
| Rosalind k=4 | `-C` | `-l 6 -h 9 -i 4` | 2 116 6 12 10 2 |
| Rosalind k=4 | plain | `-l 4 -h 4` | 3 192 4 6 5 11 |
| BA1B k=4 | `-C` | `-f -h 5` | 0 0 1 16 2 2 3 1 4 1 5 0 6 0 |
| BA1B k=4 | `-C` | `-f -l 2 -h 6` | 1 16 2 2 3 1 4 1 5 0 6 0 7 0 |
| A^31… k=3 | `-C` | (defaults) | 1 1 2 2 6 1 8 1 30 1 |
| A^31… k=3 | plain | `-h 5` | 1 6 3 2 4 2 6 1 (count 29 → cap bin 6) |
| A^31… k=3 | plain | `-f -l 3 -h 8 -i 2` | 1 6 3 4 5 0 7 0 9 1 |

Defaults pool multiplicities > 10000 into bin 10001 (`A^10010`, k=1 → `10001 1`); `--full` with defaults lists the 10002 bins 0…10001. `GetKmerSpectrum` is unchanged (no cap).

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
7. Marçais, G., Kingsford, C. (2011). Jellyfish — `sub_commands/histo_main.cc`. https://github.com/gmarcais/Jellyfish
8. scikit-bio `Sequence.kmer_frequencies`; SciPy `scipy.stats.entropy` (reference implementations used for cross-check).
