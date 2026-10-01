# Test Specification: KMER-FREQ-001

**Test Unit ID:** KMER-FREQ-001
**Area:** K-mer Analysis
**Title:** K-mer Frequency Analysis
**Created:** 2026-01-23
**Status:** Complete

---

## Canonical Methods

| Method | Class | Type | Description |
|--------|-------|------|-------------|
| `GetKmerSpectrum(string, k)` | KmerAnalyzer | Spectrum | Frequency-of-frequency distribution |
| `GetKmerFrequencies(string, k)` | KmerAnalyzer | Normalized | Frequency (0.0–1.0) per k-mer |
| `CalculateKmerEntropy(string, k)` | KmerAnalyzer | Entropy | Shannon entropy in bits |

---

## Evidence Summary

| Source | Type | Key Contributions |
|--------|------|-------------------|
| Wikipedia (K-mer) | Primary | K-mer spectrum definition, frequency distribution, genomic signatures |
| Wikipedia (Entropy) | Primary | Shannon entropy formula H = -Σ p log₂(p), max entropy = log₂(n) |
| Shannon (1948) | Primary | Original entropy definition, maximum when equiprobable |
| Rosalind KMER | Primary | K-mer composition problem with sample dataset |
| Chor et al. (2009) Genome Biol 10:R108 | Primary | k-mer spectrum = number of distinct k-mers per abundance |
| Jellyfish `histo_main.cc` | Reference impl. | `++histo[count]` per distinct k-mer; sparse output |
| scikit-bio `Sequence.kmer_frequencies`, `scipy.stats.entropy` | Reference impl. | relative freq = count/(L−k+1); H in bits with base=2 |

---

## Must Tests (M) - Evidence-Backed

### M1: GetKmerFrequencies - Frequency Sum Invariant
**Source:** Mathematical definition of probability distribution
**Test:** Sum of all k-mer frequencies equals 1.0 (within tolerance)
**Cases:**
- Standard sequence with mixed k-mers
- Homopolymer sequence
- Various k values

### M2: GetKmerFrequencies - Edge Cases
**Source:** Wikipedia pseudocode, implementation contract
**Tests:**
- Empty sequence → empty dictionary
- k > sequence length → empty dictionary

### M3: GetKmerFrequencies - Calculation Correctness
**Source:** Mathematical definition
**Test:** For sequence "AAA" with k=2: {"AA": 1.0} (only k-mer appears with 100% frequency)

### M4: GetKmerSpectrum - Spectrum Correctness
**Source:** Wikipedia K-mer: "k-mer spectrum shows the multiplicity of each k-mer"
**Test:** For "ACGTACGT" with k=4:
- ACGT appears twice, 3 others appear once
- Spectrum: {1: 3, 2: 1}

### M5: GetKmerSpectrum - Spectrum Total Invariant
**Source:** Mathematical definition
**Test:** Sum of (multiplicity × count) over spectrum entries = L - k + 1

### M6: CalculateKmerEntropy - Zero Entropy (Homopolymer)
**Source:** Shannon (1948): "When entropy is zero, there is no uncertainty"
**Test:** Homopolymer sequence (e.g., "AAAA", k=2) → entropy = 0.0

### M7: CalculateKmerEntropy - Maximum Entropy (Uniform Distribution)
**Source:** Shannon (1948): "Maximum uncertainty when all outcomes are equally likely"
**Test:** For "ACGT" with k=1 (all 4 bases once) → entropy = log₂(4) = 2.0 bits

### M8: CalculateKmerEntropy - Edge Cases
**Source:** Implementation contract
**Tests:**
- Empty sequence → 0.0
- k > sequence length → 0.0 (no k-mers to measure)

### M9: CalculateKmerEntropy - Bounds Invariant
**Source:** Shannon (1948): 0 ≤ H ≤ log₂(n)
**Test:** For any sequence, 0 ≤ entropy ≤ log₂(unique_kmer_count)

---

## Should Tests (S) - Additional Coverage

### S1: GetKmerFrequencies - Individual Frequency Correctness
**Test:** Verify specific k-mer frequencies against manual calculation

### S2: GetKmerSpectrum - Multiple Multiplicities
**Test:** Sequence with k-mers appearing at different frequencies produces correct spectrum

### S3: CalculateKmerEntropy - Intermediate Entropy
**Test:** Sequence with non-uniform k-mer distribution produces entropy between 0 and max

### S4: Case Insensitivity
**Test:** Mixed case sequences produce same results as uppercase

---

## Could Tests (C) - Extended Coverage

### C1: Performance with Larger Sequences
**Test:** Methods complete in reasonable time for longer sequences

### C2: Various K Values
**Test:** Methods work correctly across different k values (1, 2, 3, 4)

---

## Open Questions / Decisions

None. All behavior is well-defined by mathematical definitions and sources.

---

## Coverage Classification

**Test Location:** `KmerAnalyzer_Frequency_Tests.cs` (30 tests)

| Spec | Test | Classification | Notes |
|------|------|----------------|-------|
| M1 | `GetKmerFrequencies_StandardSequence_SumsToOne` | ✅ Covered | Exact sum = 1.0 |
| M1 | `GetKmerFrequencies_Homopolymer_SumsToOne` | ✅ Covered | Homopolymer variant |
| M1 | `GetKmerFrequencies_VariousK_SumsToOne` (×4) | ✅ Covered | k = 1..4 |
| M2 | `GetKmerFrequencies_EmptySequence_ReturnsEmptyDictionary` | ✅ Covered | Edge: L = 0 |
| M2 | `GetKmerFrequencies_KGreaterThanLength_ReturnsEmptyDictionary` | ✅ Covered | Edge: k > L |
| M3 | `GetKmerFrequencies_SingleKmerType_HasFrequencyOne` | ✅ Covered | count = 1, freq = 1.0 |
| S1 | `GetKmerFrequencies_MixedSequence_CorrectFrequencies` | ✅ Covered | Exact values: AA=0.4, AC/CG/GT=0.2 |
| M4 | `GetKmerSpectrum_StandardSequence_ReturnsCorrectSpectrum` | ✅ Covered | Exact: {1:3, 2:1} |
| M4 | `GetKmerSpectrum_Homopolymer_SingleEntryWithHighCount` | ✅ Covered | Exact: {4:1} |
| S2 | `GetKmerSpectrum_MultipleMultiplicities_AllCaptured` | ✅ Covered | **Strengthened**: exact values {1:1, 2:3} |
| M5 | `GetKmerSpectrum_TotalInvariant` (×4) | ✅ Covered | Σ(mult×count) = L−k+1 |
| — | `GetKmerSpectrum_EmptySequence_ReturnsEmptyDictionary` | ✅ Covered | Edge: L = 0 |
| — | `GetKmerSpectrum_KGreaterThanLength_ReturnsEmptyDictionary` | ✅ Covered | Edge: k > L |
| M6 | `CalculateKmerEntropy_Homopolymer_ReturnsZero` | ✅ Covered | H = 0.0 |
| M6 | `CalculateKmerEntropy_SingleKmer_ReturnsZero` | ✅ Covered | H = 0.0, k = L |
| M7 | `CalculateKmerEntropy_UniformDistribution_ReturnsMaxEntropy` | ✅ Covered | H = log₂(4) = 2.0 bits |
| M7 | `CalculateKmerEntropy_AllDistinctKmers_ReturnsLogOfCount` | ✅ Covered | H = log₂(3) |
| M8 | `CalculateKmerEntropy_EmptySequence_ReturnsZero` | ✅ Covered | Edge: L = 0 |
| M8 | `CalculateKmerEntropy_KGreaterThanLength_ReturnsZero` | ✅ Covered | Edge: k > L |
| M9 | `CalculateKmerEntropy_BoundsInvariant` (×4) | ✅ Covered | 0 ≤ H ≤ log₂(n) |
| S3 | `CalculateKmerEntropy_NonUniformDistribution_ExactValue` | ✅ Covered | **Strengthened**: H = log₂(5)−0.4 ≈ 1.9219 |
| S4 | `KmerFrequencyMethods_MixedCase_SameAsUppercase` | ✅ Covered | **Strengthened**: key-value pair comparison |
| C2 | `KmerFrequencyMethods_VariousKValues_InvariantsHold` | ✅ Covered | **Strengthened**: 3 invariants per k |

### Audit Log

#### Strengthening Audit (2026-03-06)

4 weak tests strengthened with exact source-backed assertions:
- **S2**: `ContainsKey` → exact `spectrum[1]=1, spectrum[2]=3`
- **S3**: range check `0 < H < max` → exact `H = log₂(5) − 0.4` (cross-verifies S1 frequencies)
- **S4**: count-only comparison → full key-value pair matching across all three methods
- **C2**: permissive `> 0` checks → three source-backed invariants (M1 + M5 + M9) per k value

---

#### Review 2026-09 (campaign B06)

Stage A PASS-with-notes, Stage B PASS. Added reference cross-check tests (region
"Reference cross-check (review 2026-09)"):
- `CalculateKmerEntropy_MatchesScipyEntropyBase2` (×5) — `scipy.stats.entropy(counts, base=2)`
- `GetKmerFrequencies_MatchesScikitBioRelativeFrequencies` — skbio `kmer_frequencies(k, relative=True)`, denominator L−k+1
- `GetKmerSpectrum_MatchesJellyfishHistoSemantics` — count-of-counts, sparse (Jellyfish `histo` without `--full`; Chor et al. 2009)
- `CalculateKmerEntropy_AgreesWithSequenceComplexityKmerEntropy` (×3) — same quantity as SEQ-COMPLEX-KMER-001
- `FrequencyMethods_NonPositiveK_NonEmptySequence_Throws`

Notes: k ≤ 0 throws only for non-empty input (empty/null short-circuits to empty/0 in `CountKmers`);
spectrum has no upper cap bin (Jellyfish `-h` default 10000 bins higher counts together) and is
single-strand (Jellyfish `-C` counts canonical k-mers); per-window cost is O(k) (substring hashing),
so total is O(n·k), not O(n).

#### Audit round 1, WP3 (B06, F12) — Jellyfish `histo` options

New method `GetKmerHistogram(sequence, k, KmerCountingOptions = default, low = 1, high = 10000, increment = 1, full = false)`
(+ `GetKmerHistogram(IEnumerable<int> counts, …)`), a verbatim port of Jellyfish `sub_commands/histo_main.cc`.
`GetKmerSpectrum` is unchanged (sparse, no cap). Tests (`KmerAnalyzer_HistogramClumpWindowsFilters_Tests`):

| ID | Test | Evidence |
|---|---|---|
| H1 | `GetKmerHistogram_MatchesJellyfishHisto` ×72 | Jellyfish 2.3.1 binary: 3 inputs × {plain, `-C`} × 12 option sets (defaults, `-h 5`, `-i 2`, `-l 3 -h 8 -i 2`, `-l 2 -h 6`, `-l 1 -h 4 -i 3`, `-l 4 -h 4`, `-f -h 5`, `-f -l 3 -h 8 -i 2`, `-f -l 2 -h 6`, `-l 5 -h 7`, `-l 6 -h 9 -i 4`) |
| H2 | `GetKmerHistogram_CountTableOverload_EqualsSequenceOverload` | table overload = sequence overload; Rosalind `-C -l 3 -h 8 -i 2` → 1 57 3 52 5 12 7 4 9 5 |
| H3 | `GetKmerHistogram_DefaultsWithoutCap_EqualSortedSpectrum` | high above all counts → non-zero spectrum bins |
| H4 | `GetKmerHistogram_DefaultHigh_PoolsCountsAboveTenThousandInCapBin10001` | `histo_main.cc`: ceil = high + inc; A^10010 k=1 → `10001 1` |
| H5 | `GetKmerHistogram_FullDefaults_Has10002BinsFromZero` | nb_buckets = (ceil + inc − base)/inc = 10002, labels 0…10001 |
| H6 | `GetKmerHistogram_FrequenciesSumToDistinctKmers` | every distinct k-mer lands in exactly one bucket |
| H7 | `GetKmerHistogram_EmptyInput_NoRows_FullStillListsBins` | `-f` on an empty table lists the zero buckets |
| H8 | `GetKmerHistogram_InvalidParameters_Throw` | high < low (Jellyfish error), inc < 1, low < 0, overflow, oversize full, null/negative counts |
| H9 | `GetKmerSpectrum_Unchanged_NoCapBin` | backward compatibility |

---

## Deviations and Assumptions

None. Implementation matches external source definitions exactly.
