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

## Audit Log

### Initial Audit (2026-01-23)

**Existing Tests Location:** `KmerAnalyzerTests.cs`

| Test | Status | Classification |
|------|--------|----------------|
| `GetKmerSpectrum_ReturnsFrequencyDistribution` | Exists | Covered (basic) |
| `GetKmerFrequencies_SumsToOne` | Exists | Covered (invariant) |
| `GetKmerFrequencies_CalculatesCorrectly` | Exists | Covered |
| `CalculateKmerEntropy_UniformDistribution_HighEntropy` | Exists | Covered |
| `CalculateKmerEntropy_SingleRepeated_ZeroEntropy` | Exists | Covered |
| `CalculateKmerEntropy_EmptySequence_ReturnsZero` | Exists | Covered |

**Gaps Identified:**
- Missing edge case tests for GetKmerSpectrum
- Missing spectrum total invariant test
- Missing entropy bounds invariant test
- Missing spectrum edge case tests
- Tests need consolidation into dedicated test file

**Consolidation Plan:**
1. Create new canonical file: `KmerAnalyzer_Frequency_Tests.cs`
2. Move relevant tests from `KmerAnalyzerTests.cs` (remove from there)
3. Add missing Must tests
4. Keep only KMER-FIND-001 related tests in `KmerAnalyzerTests.cs` for future processing

---

## ASSUMPTIONS

None. All test rationales are traceable to sources.
