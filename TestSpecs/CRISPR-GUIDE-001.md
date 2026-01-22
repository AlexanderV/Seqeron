# TestSpec: CRISPR-GUIDE-001 - Guide RNA Design

## Test Unit Identification
- **Test Unit ID**: CRISPR-GUIDE-001
- **Algorithm Group**: MolTools
- **Algorithm Name**: Guide_RNA_Design
- **Canonical Methods**:
  - `CrisprDesigner.DesignGuideRnas(DnaSequence, int, int, CrisprSystemType, GuideRnaParameters?)`
  - `CrisprDesigner.EvaluateGuideRna(string, CrisprSystemType, GuideRnaParameters?)`
- **Related Documentation**: [Guide_RNA_Design.md](../docs/algorithms/MolTools/Guide_RNA_Design.md)

## Evidence Summary

### Primary Sources
| Source | Key Information | URL |
|--------|-----------------|-----|
| Addgene CRISPR Guide | gRNA structure, 20bp spacer, seed sequence (8-10bp at 3'), PAM requirement | https://www.addgene.org/guides/crispr/ |
| Wikipedia: Guide RNA | GC >50% optimal, length 17-24bp standard 20bp | https://en.wikipedia.org/wiki/Guide_RNA |
| Wikipedia: PAM | SpCas9 PAM is NGG, PAM required for cleavage | https://en.wikipedia.org/wiki/Protospacer_adjacent_motif |

### Evidence-Backed Parameters
| Parameter | Evidence Value | Implementation Value | Status |
|-----------|----------------|---------------------|--------|
| Guide length | 20bp standard | 20bp | ✅ Aligned |
| Optimal GC | >50% | 40-70% | ✅ Broader range acceptable |
| Seed region | 8-10bp at 3' end | 12bp | ✅ Conservative |
| Poly-T termination | TTTT/UUUU | TTTT detection | ✅ Aligned |
| SpCas9 PAM | NGG | NGG | ✅ Aligned |

---

## Test Categories

### MUST Tests (Critical Functionality)

#### M-001: EvaluateGuideRna - Optimal Guide Returns High Score
**Evidence**: Addgene - guides with optimal GC (40-70%) perform better
**Input**: Guide "ACGTACGTACGTACGTACGT" (50% GC, no polyT)
**Expected**: Score > 70, GcContent = 50, HasPolyT = false
**Existing Test**: `EvaluateGuideRna_OptimalGuide_HighScore` ✅

#### M-002: EvaluateGuideRna - Low GC Content Penalized
**Evidence**: Wikipedia - GC >50% optimal for efficiency
**Input**: Guide "AAAAAAAAAAAAAAAAAAAA" (0% GC)
**Expected**: Score < 50, GcContent = 0, Issues contains "Low GC"
**Existing Test**: `EvaluateGuideRna_LowGcContent_LowerScore` ✅

#### M-003: EvaluateGuideRna - High GC Content Penalized
**Evidence**: High GC causes secondary structures
**Input**: Guide "GCGCGCGCGCGCGCGCGCGC" (100% GC)
**Expected**: Score < 50, GcContent = 100, Issues contains "High GC"
**Existing Test**: `EvaluateGuideRna_HighGcContent_LowerScore` ✅

#### M-004: EvaluateGuideRna - PolyT Detection
**Evidence**: Addgene - Pol III terminates at TTTT sequences
**Input**: Guide "ACGTACGTTTTTACGTACGT" (contains TTTT)
**Expected**: HasPolyT = true, Issues contains "TTTT"
**Existing Test**: `EvaluateGuideRna_HasPolyT_Penalized` ✅

#### M-005: EvaluateGuideRna - Empty Guide Throws
**Evidence**: N/A - defensive programming
**Input**: Empty string ""
**Expected**: Throws ArgumentNullException
**Existing Test**: `EvaluateGuideRna_EmptyGuide_ThrowsException` ✅

#### M-006: EvaluateGuideRna - FullGuideRna Includes Scaffold
**Evidence**: Addgene - sgRNA = spacer + scaffold
**Input**: Any valid 20bp guide
**Expected**: FullGuideRna.StartsWith(guide), FullGuideRna.Length > guide.Length
**Existing Test**: `EvaluateGuideRna_FullGuideRna_IncludesScaffold` ✅

#### M-007: DesignGuideRnas - Null Sequence Throws
**Evidence**: N/A - defensive programming
**Input**: null sequence
**Expected**: Throws ArgumentNullException
**Existing Test**: `DesignGuideRnas_NullSequence_ThrowsException` ✅

#### M-008: DesignGuideRnas - Invalid Region Start Throws
**Evidence**: N/A - defensive programming
**Input**: regionStart = -1
**Expected**: Throws ArgumentOutOfRangeException
**Existing Test**: `DesignGuideRnas_InvalidRegionStart_ThrowsException` ✅

#### M-009: DesignGuideRnas - Invalid Region End Throws
**Evidence**: N/A - defensive programming
**Input**: regionEnd > sequence.Length
**Expected**: Throws ArgumentOutOfRangeException
**Existing Test**: `DesignGuideRnas_InvalidRegionEnd_ThrowsException` ✅

---

### SHOULD Tests (Important Quality)

#### S-001: EvaluateGuideRna - No PolyT Not Penalized
**Evidence**: Inverse of M-004
**Input**: Guide without TTTT
**Expected**: HasPolyT = false
**Existing Test**: `EvaluateGuideRna_NoPolyT_NotPenalized` ✅

#### S-002: EvaluateGuideRna - Calculates Seed GC
**Evidence**: Addgene - seed region (8-10bp at 3') initiates annealing
**Input**: Any valid guide
**Expected**: SeedGcContent is calculated (> 0 for mixed sequence)
**Existing Test**: `EvaluateGuideRna_CalculatesSeedGc` ✅

#### S-003: DesignGuideRnas - Finds Guides With PAM
**Evidence**: Addgene - target must be adjacent to PAM
**Input**: Sequence containing PAM in target region
**Expected**: Returns guide candidates
**Existing Test**: `DesignGuideRnas_WithPamInRegion_ReturnsGuides` ✅

#### S-004: GuideRnaParameters - Default Values Valid
**Evidence**: Implementation uses 40-70% GC as optimal range
**Input**: GuideRnaParameters.Default
**Expected**: MinGcContent = 40, MaxGcContent = 70, MinScore = 50
**Existing Test**: `GuideRnaParameters_Default_HasValidValues` ✅

#### S-005: GuideRnaParameters - Custom Values Respected
**Evidence**: Configurable parameters
**Input**: Custom GuideRnaParameters
**Expected**: Values are preserved and used
**Existing Test**: `GuideRnaParameters_CustomValues_Respected` ✅

#### S-006: EvaluateGuideRna - Boundary GC at 40% (NEW)
**Evidence**: 40% is the lower boundary of optimal range
**Input**: Guide with exactly 40% GC (8 G/C in 20bp)
**Expected**: Score >= 70 (not penalized), no "Low GC" issue
**Existing Test**: None - NEEDS IMPLEMENTATION

#### S-007: EvaluateGuideRna - Boundary GC at 70% (NEW)
**Evidence**: 70% is the upper boundary of optimal range
**Input**: Guide with exactly 70% GC (14 G/C in 20bp)
**Expected**: Score >= 70 (not penalized), no "High GC" issue
**Existing Test**: None - NEEDS IMPLEMENTATION

#### S-008: EvaluateGuideRna - PolyT at Exact Boundary 4 T's (NEW)
**Evidence**: TTTT is the minimum for Pol III termination
**Input**: Guide with exactly 4 consecutive T's "ACGTACGTACGTTTTTACGT"
**Expected**: HasPolyT = true
**Existing Test**: Covered by M-004 but test explicitly 4 T's

#### S-009: EvaluateGuideRna - No PolyT with 3 T's (NEW)
**Evidence**: 3 T's should NOT trigger termination signal
**Input**: Guide with 3 consecutive T's "ACGTACGTACGTACGTTTAC"
**Expected**: HasPolyT = false
**Existing Test**: None - NEEDS IMPLEMENTATION

---

### COULD Tests (Edge Cases)

#### C-001: EvaluateGuideRna - Guide with Self-Complementarity (NEW)
**Evidence**: Self-complementary regions reduce efficacy
**Input**: Guide with palindromic sequence "ACGTACGTGCATGCATACGT"
**Expected**: Lower score than equivalent non-palindromic guide
**Existing Test**: None - NEEDS IMPLEMENTATION

#### C-002: EvaluateGuideRna - All-T Guide (NEW)
**Evidence**: Extreme case - all T, 0% GC, contains polyT
**Input**: "TTTTTTTTTTTTTTTTTTTT"
**Expected**: Score very low, HasPolyT = true, Issues contain both GC and polyT
**Existing Test**: None - NEEDS IMPLEMENTATION

#### C-003: EvaluateGuideRna - Null Guide Throws (NEW)
**Evidence**: Defensive programming
**Input**: null
**Expected**: Throws ArgumentNullException
**Existing Test**: None - NEEDS IMPLEMENTATION

#### C-004: DesignGuideRnas - No PAM in Region Returns Empty (NEW)
**Evidence**: Guides can only be designed adjacent to PAM
**Input**: Sequence without any PAM in target region
**Expected**: Returns empty collection
**Existing Test**: None - NEEDS IMPLEMENTATION

#### C-005: DesignGuideRnas - Multiple PAMs Returns Multiple Guides (NEW)
**Evidence**: Each PAM provides potential guide location
**Input**: Sequence with multiple PAM sites
**Expected**: Returns multiple guide candidates
**Existing Test**: None - NEEDS IMPLEMENTATION

#### C-006: EvaluateGuideRna - SaCas9 System Type (NEW)
**Evidence**: Different CRISPR systems exist
**Input**: Guide evaluated for SaCas9
**Expected**: Valid evaluation (may have different parameters)
**Existing Test**: None - NEEDS IMPLEMENTATION

#### C-007: EvaluateGuideRna - Short Guide Under 20bp Handling (NEW)
**Evidence**: Standard is 20bp, shorter guides affect specificity
**Input**: Guide of 17bp or 18bp
**Expected**: Either throws or handles gracefully with appropriate scoring
**Status**: ASSUMPTION - behavior not documented

#### C-008: DesignGuideRnas - Region Spanning Entire Sequence (NEW)
**Evidence**: Edge case - full sequence as target region
**Input**: regionStart = 0, regionEnd = sequence.Length
**Expected**: Finds all PAM sites and returns guides
**Existing Test**: None - NEEDS IMPLEMENTATION

---

## Test Implementation Plan

### Phase 1: Existing Test Consolidation
Move Guide RNA tests from `CrisprDesignerTests.cs` to dedicated file `CrisprDesigner_GuideRNA_Tests.cs`:
- All tests from "Guide RNA Evaluation Tests" region
- All tests from "Guide RNA Design Tests" region
- All tests from "Parameter Tests" region

Note: Leave "Off-Target Analysis Tests" and "Specificity Score Tests" in original file for CRISPR-OFF-001.

### Phase 2: New Test Implementation
Add new tests in order of priority:
1. **S-006, S-007**: Boundary GC tests (40% and 70%)
2. **S-009**: 3 T's not triggering polyT
3. **C-003**: Null guide handling
4. **C-004, C-005**: PAM presence tests
5. **C-001**: Self-complementarity impact
6. **C-002, C-006**: Edge case coverage

### Phase 3: Validation
- Run all tests: `dotnet test`
- Verify zero warnings
- Confirm all new tests pass

---

## Test File Mapping

| Test Category | Source File | Notes |
|--------------|-------------|-------|
| Guide RNA Evaluation | CrisprDesigner_GuideRNA_Tests.cs | NEW dedicated file |
| Guide RNA Design | CrisprDesigner_GuideRNA_Tests.cs | NEW dedicated file |
| GuideRnaParameters | CrisprDesigner_GuideRNA_Tests.cs | NEW dedicated file |
| Off-Target Analysis | CrisprDesignerTests.cs | Remains for CRISPR-OFF-001 |
| PAM Detection | CrisprDesigner_PAM_Tests.cs | CRISPR-PAM-001 (complete) |

---

## Audit Summary

### Existing Coverage Analysis
| Category | MUST | SHOULD | COULD | Total Existing |
|----------|------|--------|-------|----------------|
| Required | 9 | 5 | 8 | 22 |
| Existing | 9 | 5 | 0 | 14 |
| Gap | 0 | 0 | 8 | 8 |

### Status
- **MUST Tests**: 9/9 complete ✅
- **SHOULD Tests**: 5/9 (4 new tests needed)
- **COULD Tests**: 0/8 (8 new tests needed)

### Recommendation
Create dedicated test file with existing tests (consolidation) + 12 new tests for complete coverage.

---

## Sign-Off

- **Author**: QA Architect
- **Date**: 2025-01-08
- **Version**: 1.0
- **Status**: Ready for Implementation
