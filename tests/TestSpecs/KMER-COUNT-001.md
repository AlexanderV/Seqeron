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

1. **Wikipedia - K-mer:** Definition, algorithm pseudocode, L − k + 1 formula, applications
2. **Rosalind - K-mer Composition (KMER):** 4-mer composition example with expected counts
3. **Rosalind - Clump Finding (BA1E):** K-mer clump definition for validation

## Test Categories

### MUST Tests (Evidence-Backed)

| ID | Test | Evidence |
|----|------|----------|
| M1 | Empty sequence returns empty dictionary | Wikipedia pseudocode, implementation pattern |
| M2 | k > sequence length returns empty dictionary | Wikipedia: L − k + 1 formula (negative yields 0) |
| M3 | k ≤ 0 throws ArgumentOutOfRangeException | Implementation validation |
| M4 | null sequence returns empty dictionary | Defensive programming |
| M5 | Total count invariant: sum(counts) = L − k + 1 | Wikipedia: "L − k + 1 k-mers" |
| M6 | Homopolymer: single k-mer, count = L − k + 1 | Wikipedia: all same bases example |
| M7 | Case-insensitive counting | Implementation-specific normalization |
| M8 | Distinct k-mers counted correctly | Rosalind KMER problem |
| M9 | Overlapping k-mers counted correctly | Wikipedia sliding window |
| M10 | CountKmersSpan produces same results as CountKmers | API consistency |
| M11 | CountKmersBothStrands combines forward + reverse complement | DNA double-strand property |

### SHOULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | Mixed case input normalized | Robustness for real-world data |
| S2 | Non-DNA characters handled (IUPAC N) | Genomic data often contains N |
| S3 | k = 1 counts individual nucleotides | Edge case at minimum valid k |
| S4 | k = sequence length yields single k-mer | Boundary condition |
| S5 | Large sequence performance acceptable | Practical use case |

### COULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | Cancellation token stops operation | Async API contract |
| C2 | Progress reporting works | Async API contract |

## Consolidation Plan

### Current Test Pool

| File | Tests | Status |
|------|-------|--------|
| KmerAnalyzerTests.cs | 31 tests covering multiple methods | Mixed (includes non-KMER-COUNT-001 tests) |
| PerformanceExtensionsTests.cs | 4+ tests for CountKmersSpan | Delegate smoke tests |

### Consolidation Actions

1. **Create** `KmerAnalyzer_CountKmers_Tests.cs` for KMER-COUNT-001 canonical deep tests
2. **Move** CountKmers tests from KmerAnalyzerTests.cs to new file
3. **Retain** other k-mer tests (spectrum, entropy, clumps) in KmerAnalyzerTests.cs for future test units
4. **Retain** PerformanceExtensionsTests as smoke tests for span-based API

### Existing Test Audit

| Test | Classification | Action |
|------|----------------|--------|
| `CountKmers_SimpleSequence_CountsCorrectly` | Covered | Keep, enhance |
| `CountKmers_AllSame_SingleCount` | Covered (M6) | Keep |
| `CountKmers_EmptySequence_ReturnsEmptyDictionary` | Covered (M1) | Keep |
| `CountKmers_KLargerThanSequence_ReturnsEmptyDictionary` | Covered (M2) | Keep |
| `CountKmers_InvalidK_ThrowsException` | Covered (M3) | Keep |
| `CountKmers_DnaSequence_Works` | Covered (wrapper smoke) | Keep |
| `CountKmersSpan_ReturnsCorrectCounts` | Covered (M10) | Keep in PerformanceExtensionsTests |

### Missing Tests

- M5: Total count invariant (add property-based test)
- M7: Case-insensitive counting (add explicit test)
- M11: CountKmersBothStrands validation (add explicit test)
- S3: k=1 edge case
- S4: k=sequence length edge case

## Open Questions / Decisions

None - algorithm behavior is well-defined by evidence.

## ASSUMPTIONS

None - all test rationale backed by evidence sources.
