# Test Specification: ANNOT-ORF-001

**Test Unit ID:** ANNOT-ORF-001
**Area:** Annotation
**Algorithm:** ORF Detection
**Status:** Draft
**Last Updated:** 2026-01-23

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

### GenomeAnnotatorTests.cs (Current Location)
| Test | Status | Assessment |
|------|--------|------------|
| FindOrfs_SimpleOrf_FindsIt | Covered | Weak assertions |
| FindOrfs_NoStartCodon_ReturnsEmpty | Covered | OK |
| FindOrfs_NoStopCodon_ReturnsEmpty | Covered | OK |
| FindOrfs_MultipleFrames_FindsAll | Covered | Weak - only checks frame 1 |
| FindOrfs_ReverseStrand_FinitsOrfs | Covered | Typo in name |
| FindOrfs_BelowMinLength_Excluded | Covered | OK |
| FindOrfs_EmptySequence_ReturnsEmpty | Covered | OK |
| FindOrfs_AlternativeStartCodons_Detected | Covered | OK |
| FindLongestOrfsPerFrame_ReturnsFramesDictionary | Covered | OK |
| FindLongestOrfsPerFrame_BothStrands_IncludesNegativeFrames | Covered | OK |

### TranslatorTests.cs (Wrapper)
| Test | Status | Assessment |
|------|--------|------------|
| FindOrfs_SimpleOrf_FindsIt | Duplicate | Remove or keep as smoke |
| FindOrfs_NoStartCodon_ReturnsEmpty | Duplicate | Remove or keep as smoke |
| FindOrfs_ShortOrf_FilteredByMinLength | Duplicate | Remove or keep as smoke |
| FindOrfs_RespectMinLength_FindsSmallOrfs | Duplicate | Keep as smoke |
| FindOrfs_ForwardOnly_DoesNotSearchReverseStrand | Covered | Keep as smoke |
| FindOrfs_ResultHasCorrectFrame | Covered | Keep as smoke |
| FindOrfs_NullDna_ThrowsException | Covered | Keep as smoke |

### GenomicAnalyzerTests.cs (Alternate)
| Test | Status | Assessment |
|------|--------|------------|
| FindOpenReadingFrames_SimpleOrf_FindsIt | Duplicate | Keep as smoke |
| FindOpenReadingFrames_MultipleFrames_FindsAll | Duplicate | Keep as smoke |
| FindOpenReadingFrames_NoOrf_ReturnsEmpty | Duplicate | Keep as smoke |

---

## Consolidation Plan

### Canonical Test File
**File:** `GenomeAnnotator_ORF_Tests.cs`
- Refactor existing tests from GenomeAnnotatorTests.cs
- ~~Add missing Must tests~~ ✅ Done
- Strong invariant assertions

### Wrapper Smoke Tests (Keep Minimal)
- **TranslatorTests.cs**: Keep 2-3 delegation smoke tests
- **GenomicAnalyzerTests.cs**: Keep existing 3 tests as smoke

### Tests to Remove/Refactor
- Rename typo: `FindOrfs_ReverseStrand_FinitsOrfs` → `FindOrfs_ReverseStrand_FindsOrfs`
- Strengthen weak assertions in canonical tests

---

## Open Questions

| # | Question | Decision |
|---|----------|----------|
| 1 | Behavior when ORF extends to end without stop codon? | Document current impl behavior |

---

## ASSUMPTIONS

| # | Assumption | Rationale |
|---|------------|-----------|
| A01 | N characters are skipped/ignored | No explicit source; implementation-specific |
