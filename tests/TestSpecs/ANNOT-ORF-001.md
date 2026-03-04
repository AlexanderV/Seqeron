# Test Specification: ANNOT-ORF-001

**Test Unit ID:** ANNOT-ORF-001
**Area:** Annotation
**Algorithm:** ORF Detection
**Status:** Complete
**Last Updated:** 2026-01-24

---

## Scope

### Canonical Methods (Deep Testing)
| Method | Class | Priority |
|--------|-------|----------|
| `FindOrfs(dna, minLen, bothStrands, requireStart)` | GenomeAnnotator | Canonical |
| `FindLongestOrfsPerFrame(dna, bothStrands)` | GenomeAnnotator | Per-frame variant |

### Alternate/Wrapper Methods (Smoke Testing Only)
| Method | Class | Priority |
|--------|-------|----------|
| `FindOpenReadingFrames(seq, minLen)` | GenomicAnalyzer | Alternate (smoke) |
| `FindOrfs(dna, geneticCode, minLen, bothStrands)` | Translator | Wrapper (smoke) |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia (Open reading frame) | Encyclopedia | Definition, six-frame translation, start/stop codons |
| Rosalind ORF Problem | Educational | Problem definition, sample dataset with expected outputs |
| NCBI ORF Finder | Tool Reference | Parameters (minLength, genetic code, start codon options) |
| Deonier et al. (2005) | Textbook | MinLength 150 codons typical |
| Claverie (1997) | Paper | MinLength 100 codons typical |

---

## Test Requirements

### Must Tests (Evidence-Based)

| ID | Test Case | Evidence | Expected Result |
|----|-----------|----------|-----------------|
| M01 | Simple ORF (ATG...TAA) | Wikipedia, Rosalind | ORF detected with correct start/end |
| M02 | Empty sequence | Standard edge case | Return empty collection |
| M03 | No start codon (requireStart=true) | Algorithm definition | Return empty |
| M04 | No stop codon | Algorithm definition | Return empty (or truncated per impl) |
| M05 | Alternative start codons (GTG, TTG) | Wikipedia, NCBI | ORFs detected with GTG/TTG starts |
| M06 | Minimum length filtering | NCBI, Claverie | ORFs below minLength excluded |
| M07 | Six-frame search | Wikipedia, Rosalind | ORFs found in all 6 frames |
| M08 | Reverse complement strand | Rosalind | ORFs on reverse strand detected |
| M09 | Frame 1, 2, 3 distinction | Wikipedia | Each frame reported separately |
| M10 | Multiple ORFs in sequence | General | All qualifying ORFs returned |
| M11 | Overlapping ORFs (same frame) | General | Both ORFs reported |
| M12 | Rosalind sample dataset | Rosalind | All 4 expected proteins found |
| M13 | ORF invariant: starts with start codon | Algorithm invariant | All returned ORFs start with ATG/GTG/TTG |
| M14 | ORF invariant: ends with stop codon | Algorithm invariant | All returned ORFs end with TAA/TAG/TGA |
| M15 | ORF invariant: length divisible by 3 | Algorithm invariant | All ORF DNA lengths % 3 == 0 |
| M16 | FindLongestOrfsPerFrame: returns all 6 keys | Per-frame requirement | Keys 1,2,3,-1,-2,-3 present |
| M17 | FindLongestOrfsPerFrame: correct longest | Per-frame requirement | Returns longest ORF per frame |

### Should Tests (Important Coverage)

| ID | Test Case | Rationale | Expected Result |
|----|-----------|-----------|-----------------|
| S01 | Lowercase input handling | Usability | Case-insensitive detection |
| S02 | Mixed case input | Usability | Correctly parsed |
| S03 | Very long sequence (10kb+) | Performance | Completes in reasonable time |
| S04 | Sequence with N characters | Real data | Graceful handling |
| S05 | Nested ORFs (inner start within outer) | Algorithm behavior | Both reported if qualifying |

### Could Tests (Nice to Have)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| C01 | Performance benchmark | Baseline recording |
| C02 | Alternative genetic codes | Bacterial code 11 |

---

## Audit of Existing Tests

### GenomeAnnotator_ORF_Tests.cs (Canonical — All Tests)
| Test | Spec ID | Status |
|------|---------|--------|
| FindOrfs_SimpleAtgTaaOrf_DetectsOrf | M01 | ✅ Strong |
| FindOrfs_EmptySequence_ReturnsEmpty | M02 | ✅ Strong |
| FindOrfs_NoStartCodon_RequireStart_ReturnsEmpty | M03 | ✅ Strong |
| FindOrfs_NoStopCodon_RequireStart_ReturnsEmpty | M04 | ✅ Strong |
| FindOrfs_AlternativeStartCodons_Detected (GTG, TTG) | M05 | ✅ Strong |
| FindOrfs_BelowMinLength_Excluded | M06 | ✅ Strong |
| FindOrfs_ExactlyMinLength_Included | M06b | ✅ Strong |
| FindOrfs_SixFrameSearch_FindsOrfsInMultipleFrames | M07 | ✅ Strong (single-sequence, exact frame assertions) |
| FindOrfs_ReverseStrand_FindsOrfs | M08 | ✅ Strong (exact count, protein) |
| FindOrfs_ForwardOnly_DoesNotSearchReverse | M08b | ✅ Strong |
| FindOrfs_FrameNumber_CorrectlyAssigned (0→1, 1→2, 2→3) | M09 | ✅ Strong (TestCase, exact frame + start) |
| FindOrfs_MultipleOrfs_AllReturned | M10 | ✅ Strong (exact count=2, positions) |
| FindOrfs_NestedOrfs_BothReportedIfQualifying | M11/S05 | ✅ Strong (exact count=2, positions, shared stop) |
| FindOrfs_RosalindDataset_FindsExpectedProteins | M12 | ✅ Strong |
| FindOrfs_Invariant_StartsWithStartCodon | M13 | ✅ Strong |
| FindOrfs_Invariant_EndsWithStopCodon | M14 | ✅ Strong |
| FindOrfs_Invariant_LengthDivisibleBy3 | M15 | ✅ Strong |
| FindLongestOrfsPerFrame_BothStrands_Returns6Keys | M16 | ✅ Strong |
| FindLongestOrfsPerFrame_ForwardOnly_Returns3Keys | M16b | ✅ Strong |
| FindLongestOrfsPerFrame_ReturnsLongestPerFrame | M17 | ✅ Strong |
| FindOrfs_LowercaseInput_HandledCorrectly | S01 | ✅ Strong (compared with uppercase) |
| FindOrfs_MixedCaseInput_HandledCorrectly | S02 | ✅ Strong (exact protein MKK) |
| FindOrfs_AllStopCodons_Recognized (TAA, TAG, TGA) | M14+ | ✅ Strong (supplements M14) |
| FindOrfs_VeryLongSequence_CompletesCorrectly | S03 | ✅ Strong (10kb+ sequence, exact ORF size) |
| FindOrfs_NInStartCodon_NotRecognizedAsStart | S04a | ✅ Strong |
| FindOrfs_NInStopCodon_NotRecognizedAsStop | S04b | ✅ Strong (exact count=1) |
| FindOrfs_NInCodingRegion_OrfContinues | S04c | ✅ Strong (exact count=1) |
| FindOrfs_VeryShortSequence_ReturnsEmpty | Edge | ✅ Strong |
| FindOrfs_OnlyStartCodon_ReturnsEmpty | Edge | ✅ Strong |
| FindOrfs_NullSequence_ReturnsEmpty | Edge | ✅ Strong |

### TranslatorTests.cs (Wrapper — 3 Smoke Tests)
| Test | Status |
|------|--------|
| FindOrfs_RespectMinLength_FindsSmallOrfs | ✅ Smoke (delegation) |
| FindOrfs_ForwardOnly_DoesNotSearchReverseStrand | ✅ Smoke (parameter forwarding) |
| FindOrfs_NullDna_ThrowsException | ✅ Smoke (error handling) |

### GenomicAnalyzerTests.cs (Alternate — 3 Smoke Tests)
| Test | Status |
|------|--------|
| FindOpenReadingFrames_SimpleOrf_FindsIt | ✅ Smoke |
| FindOpenReadingFrames_MultipleFrames_FindsAll | ✅ Smoke |
| FindOpenReadingFrames_NoOrf_ReturnsEmpty | ✅ Smoke |

---

## Consolidation Plan

### Completed ✅
- **Canonical Test File:** `GenomeAnnotator_ORF_Tests.cs` — 30 tests (including TestCase variants), all with strong exact assertions
- **Wrapper Smoke Tests:** TranslatorTests.cs trimmed to 3 smoke tests (removed 4 duplicates)
- **Alternate Smoke Tests:** GenomicAnalyzerTests.cs — 3 smoke tests kept as-is
- **Old tests:** ORF tests removed from GenomeAnnotatorTests.cs (consolidated into canonical file)
- **Typo fixed:** `FinitsOrfs` → `FindsOrfs`
- **All weak assertions strengthened:** `GreaterThan(0)` / `GreaterThanOrEqualTo` replaced with exact values
- **Missing test implemented:** S03 (10kb+ sequence)

---

## Deviations and Assumptions

None. All behaviors are verified against external sources.

| Item | Status | Evidence |
|------|--------|----------|
| N character handling | Resolved | NCBI C++ Toolkit `orf.cpp`: codons containing N do not match start/stop patterns. Our implementation matches this behavior — `.ToUpperInvariant()` + set lookup naturally excludes N-containing codons. |
| ORF without stop codon | Resolved | Rosalind, NCBI ORF Finder, orfipy: standard ORF detection requires a stop codon. Our implementation returns empty when no stop codon is found (with `requireStartCodon=true`). |
