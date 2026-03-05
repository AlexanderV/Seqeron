# KMER-COUNT-001: K-mer Counting Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | KMER-COUNT-001 |
| **Area** | K-mer Analysis |
| **Canonical Method** | `KmerAnalyzer.CountKmers(string, int)` |
| **Complexity** | O(n) |
| **Invariant** | Sum of counts = n − k + 1 |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `CountKmers(string, int)` | KmerAnalyzer | Canonical | Deep |
| `CountKmersSpan(ReadOnlySpan<char>, int)` | SequenceExtensions | Span variant | Deep |
| `CountKmersBothStrands(DnaSequence, int)` | KmerAnalyzer | Both strands | Deep |
| `CountKmers(DnaSequence, int)` | KmerAnalyzer | Wrapper | Smoke |
| `CountKmers(string, int, CancellationToken, IProgress<double>)` | KmerAnalyzer | Async delegate | Smoke |

## Evidence Sources

1. **Wikipedia — K-mer:** Definition, pseudocode, L − k + 1 formula, k-mer tables, 4^k bound
   - URL: https://en.wikipedia.org/wiki/K-mer
2. **Rosalind — K-mer Composition (KMER):** 4-mer composition sample dataset (415 bp, 209 unique 4-mers, sum = 412)
   - URL: https://rosalind.info/problems/kmer/
3. **Rosalind — Clump Finding (BA1E):** K-mer clump definition
   - URL: https://rosalind.info/problems/ba1e/

## Test Categories

### MUST Tests (Evidence-Backed)

| ID | Test | Evidence |
|----|------|----------|
| M1 | Empty sequence returns empty dictionary | Wikipedia pseudocode: loop 0 to L−k+1 yields nothing when L=0 |
| M2 | k > sequence length returns empty dictionary | Wikipedia: L − k + 1 becomes ≤ 0 |
| M3 | k ≤ 0 throws ArgumentOutOfRangeException (all APIs) | k must be positive for valid k-mer definition |
| M4 | null sequence returns empty dictionary | Defensive programming |
| M5 | Total count invariant: sum(counts) = L − k + 1 | Wikipedia: "a sequence of length L will have L − k + 1 k-mers" |
| M6 | Homopolymer: single k-mer, count = L − k + 1 | Wikipedia: all same bases example |
| M7 | Case-insensitive counting (all APIs) | Wikipedia/algorithm norm: k-mers case-insensitive |
| M8 | Distinct k-mers counted correctly | Rosalind KMER problem |
| M9 | Overlapping k-mers counted correctly | Wikipedia sliding window pseudocode |
| M10 | CountKmersSpan produces same results as CountKmers (including case) | API consistency |
| M11 | CountKmersBothStrands combines forward + reverse complement | DNA double-strand property |
| M12 | Wikipedia example: "ATGG" → ATG, TGG | Wikipedia lead diagram caption |
| M13 | Wikipedia table: "GTAGAGCTGT" exact k-mers for k=2 and k=4 | Wikipedia Introduction table |
| M14 | Rosalind KMER sample: specific 4-mer counts + exact unique count (209) | Rosalind sample output |
| M15 | CancellationToken overload normalizes case | API consistency with canonical method |

### SHOULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | Mixed case input normalized (exact counts) | Robustness for real-world data |
| S2 | Non-DNA characters handled (IUPAC N, exact counts) | Genomic data often contains N |
| S3 | k = 1 counts individual nucleotides | Edge case at minimum valid k |
| S4 | k = sequence length yields single k-mer | Boundary condition |

### COULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | Cancellation token stops operation | Async API contract |
| C2 | Progress reporting works | Async API contract |

## Test File

| File | Tests | Coverage |
|------|-------|----------|
| KmerAnalyzer_CountKmers_Tests.cs | 36 tests | All MUST/SHOULD tests covered |

## Coverage Classification

All tests use exact values derived from theory (Wikipedia/Rosalind), not from implementation output.

| Test Method | Classification | Notes |
|-------------|---------------|-------|
| CountKmers_EmptySequence_ReturnsEmptyDictionary | ✅ Covered | M1 |
| CountKmers_NullSequence_ReturnsEmptyDictionary | ✅ Covered | M4 |
| CountKmers_KLargerThanSequence_ReturnsEmptyDictionary | ✅ Covered | M2 |
| CountKmers_InvalidK_ThrowsArgumentOutOfRangeException | ✅ Covered | M3 (k=0, -1, -10) |
| CountKmers_KEqualSequenceLength_ReturnsSingleKmer | ✅ Covered | S4 |
| CountKmers_TotalCountInvariant_SumEqualsLMinusKPlusOne | ✅ Covered | M5 (4 cases) |
| CountKmers_TotalCountInvariant_HoldsForAllValidK | ✅ Covered | M5 (k=1..12) |
| CountKmers_SimpleSequence_CountsDistinctKmersCorrectly | ✅ Covered | M8 |
| CountKmers_Homopolymer_SingleKmerWithCorrectCount | ✅ Covered | M6 (4 cases) |
| CountKmers_OverlappingKmers_AllCounted | ✅ Covered | M9 |
| CountKmers_KEqualsOne_CountsNucleotides | ✅ Covered | S3 |
| CountKmers_LowercaseSequence_NormalizedToUppercase | ✅ Covered | M7 (exact counts) |
| CountKmers_MixedCase_TreatedAsSameKmer | ✅ Covered | S1 (all k-mers verified) |
| CountKmers_WithAmbiguousBase_CountedAsIs | ✅ Covered | S2 (exact counts) |
| CountKmersSpan_ProducesSameResultAsCountKmers | ✅ Covered | M10 |
| CountKmersSpan_EmptySpan_ReturnsEmptyDictionary | ✅ Covered | M10 edge |
| CountKmersSpan_KLargerThanSpan_ReturnsEmptyDictionary | ✅ Covered | M10 edge |
| CountKmersSpan_InvalidK_ThrowsArgumentOutOfRangeException | ✅ Covered | M3 (Span API, k=0, -1) |
| CountKmersSpan_KEqualToLength_ReturnsSingleKmer | ✅ Covered | M10 boundary |
| CountKmersBothStrands_CombinesForwardAndReverseComplement | ✅ Covered | M11 palindromic |
| CountKmersBothStrands_NonPalindromicSequence_AddsNewKmers | ✅ Covered | M11 non-palindromic |
| CountKmersBothStrands_TotalCountInvariant | ✅ Covered | M11 invariant |
| CountKmers_DnaSequence_DelegatesToStringVersion | ✅ Covered | Wrapper smoke |
| CountKmers_WikipediaExample_ATGG_TwoThreeMers | ✅ Covered | M12 |
| CountKmers_WikipediaTable_GTAGAGCTGT_TwoMers | ✅ Covered | M13 k=2 (all 7 k-mers) |
| CountKmers_WikipediaTable_GTAGAGCTGT_TotalKmersPerK | ✅ Covered | M13 k=1..10 |
| CountKmers_WikipediaTable_GTAGAGCTGT_FourMers | ✅ Covered | M13 k=4 (all 7 k-mers) |
| CountKmers_RosalindKmerSample_SpecificCounts | ✅ Covered | M14 (7 spot checks + total) |
| CountKmers_RosalindKmerSample_ExactUniqueCount | ✅ Covered | M14 (exact = 209) |
| CountKmersSpan_LowercaseInput_NormalizesToUppercase | ✅ Covered | M7/M10 regression |
| CountKmersSpan_MixedCase_MatchesCountKmers | ✅ Covered | M10 case regression |
| CountKmers_CancellationOverload_NormalizesCase | ✅ Covered | M15 |

## Deviations and Assumptions

None — all behavior verified against external evidence sources.
