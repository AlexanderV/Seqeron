# Test Specification: REP-DIRECT-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-DIRECT-001 |
| **Area** | Repeats |
| **Title** | Direct Repeat Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-03-02 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindDirectRepeats(DnaSequence, minLength, maxLength, minSpacing)` | RepeatFinder | Canonical | Deep testing |
| `FindDirectRepeats(string, minLength, maxLength, minSpacing)` | RepeatFinder | Overload | Smoke testing |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - Direct repeat](https://en.wikipedia.org/wiki/Direct_repeat) | Definition | Sequence repeated with same directionality; may have intervening nucleotides |
| [Wikipedia - Repeated sequence (DNA)](https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)) | Context | Direct vs inverted repeats; types: tandem, interspersed, flanking |
| Ussery et al. (2009) | Technical | Computing for Comparative Microbial Genomics, Springer, Chapter 8 |
| Richard (2021) PMC8145212 | Clinical | Trinucleotide repeat expansions and mismatch repair |

---

## Test Categories

### MUST Tests (Required for DoD)

All MUST tests are justified by evidence or explicitly marked.

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| M1 | SimpleDirectRepeat_FindsRepeat | Core algorithm - detects identical sequences at two positions | Wikipedia - Direct repeat |
| M2 | AdjacentRepeats_WithZeroSpacing_Found | Adjacent repeats (minSpacing=0) should be detected | Wikipedia - tandem direct repeats |
| M3 | NoRepeats_ReturnsEmpty | Sequence without repeated patterns | Standard edge case |
| M4 | EmptySequence_ReturnsEmpty | Boundary - empty input | Standard boundary |
| M5 | NullSequence_ThrowsArgumentNullException | Parameter validation | Implementation contract |
| M6 | MinLengthTooSmall_ThrowsException | minLength < 2 is invalid | Implementation contract |
| M7 | MaxLengthLessThanMinLength_ThrowsException | Invalid parameter combination | Implementation contract |
| M8 | SpacingCalculation_Correct | Spacing = SecondPosition - FirstPosition - Length | Invariant |
| M9 | FirstPosition_LessThanSecondPosition | FirstPosition < SecondPosition always | Invariant |
| M10 | RepeatSequence_MatchesActualSequence | RepeatSequence equals substring at FirstPosition | Invariant |
| M11 | MinLength_RespectsThreshold | Only repeats ≥ minLength returned | Algorithm specification |
| M12 | MaxLength_RespectsThreshold | Only repeats ≤ maxLength returned | Algorithm specification |
| M13 | MinSpacing_RespectsThreshold | Only repeats with spacing ≥ minSpacing returned | Algorithm specification |
| M14 | SequenceTooShort_ReturnsEmpty | Sequence shorter than 2×minLength | Boundary condition |

### SHOULD Tests (Important but not blocking)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| S1 | MultipleDirectRepeats_FindsAllPairs | Three occurrences produce all three pairwise combinations | Wikipedia - "nucleotide sequences present in multiple copies" |
| S2 | StringOverload_MatchesDnaSequenceOverload | API consistency | Implementation contract |
| S3 | CaseInsensitivity_HandledCorrectly | Lowercase input processed | Implementation robustness |
| S4 | LongSpacing_Detected | Repeats with large intervening region | Wikipedia - interspersed repeats |
| S5 | BiologicalRepeat_TrinucleotideCAG | Test with disease-relevant repeat | Richard (2021) |

### COULD Tests (Nice to have)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| C1 | OverlappingPatterns_HandledCorrectly | Multiple patterns at same position | Edge case |
| C2 | LargeSequence_Performance | O(n²) completes in bounded time for 1000bp | DoD requirement |

---

## Test Audit

**Canonical file:** `RepeatFinder_DirectRepeat_Tests.cs` — 21 tests (14 MUST, 5 SHOULD, 2 COULD).

---

## Deviations and Assumptions

None. All behaviors and test data verified against external sources (Wikipedia: Direct repeat, Repeated sequence (DNA); Ussery et al. 2009; Richard 2021 PMC8145212), accessed 2026-03-02.
