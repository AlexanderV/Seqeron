# TestSpec: TRANS-PROT-001 - Protein Translation

**Test Unit ID:** TRANS-PROT-001  
**Area:** Translation  
**Canonical Class:** `Translator`  
**Date:** 2026-02-04

---

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `Translate(DnaSequence, geneticCode, frame, toFirstStop)` | Translator | Canonical | Deep |
| `Translate(RnaSequence, geneticCode, frame, toFirstStop)` | Translator | Canonical | Deep |
| `Translate(string, geneticCode, frame, toFirstStop)` | Translator | Canonical | Deep |
| `TranslateSixFrames(DnaSequence, geneticCode)` | Translator | Canonical | Deep |
| `FindOrfs(DnaSequence, geneticCode, minLength, searchBothStrands)` | Translator | Canonical | Deep |

---

## Test Classification

### Must Tests (Evidence-backed)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M01 | Translate_SingleCodon_ReturnsCorrectAminoAcid | Basic translation | Wikipedia: Translation |
| M02 | Translate_MultipleCodens_ConcatenatesAminoAcids | Sequential reading | Wikipedia: Translation |
| M03 | Translate_ToFirstStop_TerminatesAtStopCodon | Stop codon behavior | Wikipedia: Translation |
| M04 | Translate_Frame1_ShiftsReadingByOne | Frame offset mechanism | Wikipedia: Reading frame |
| M05 | Translate_Frame2_ShiftsReadingByTwo | Frame offset mechanism | Wikipedia: Reading frame |
| M06 | Translate_InvalidFrame_ThrowsException | Frame validation | Implementation spec |
| M07 | Translate_EmptySequence_ReturnsEmpty | Edge case | Implementation spec |
| M08 | Translate_NullDna_ThrowsArgumentNullException | Null handling | Implementation spec |
| M09 | Translate_NullRna_ThrowsArgumentNullException | Null handling | Implementation spec |
| M10 | Translate_RnaSequence_TranslatesCorrectly | RNA input support | Wikipedia: Translation |
| M11 | Translate_DnaSequence_ConvertsTToU | DNA normalization | Wikipedia: Genetic code |
| M12 | Translate_Lowercase_CaseInsensitive | Input flexibility | Implementation spec |
| M13 | TranslateSixFrames_ReturnsSixFrames | Six-frame translation | Wikipedia: ORF |
| M14 | TranslateSixFrames_ForwardFramesMatchDirect | Consistency | Wikipedia: Reading frame |
| M15 | TranslateSixFrames_NegativeFramesUseReverseComplement | Reverse strand | Wikipedia: Reading frame |
| M16 | FindOrfs_SimpleOrf_DetectsIt | ORF detection | Wikipedia: ORF |
| M17 | FindOrfs_NoStartCodon_ReturnsEmpty | ORF requirement | Wikipedia: ORF |
| M18 | FindOrfs_MinLengthFilter_FiltersShortOrfs | Length filtering | Wikipedia: ORF |
| M19 | FindOrfs_NullInput_ThrowsException | Null handling | Implementation spec |
| M20 | Translate_AlternativeGeneticCode_UsesCorrectMappings | Alternative codes | NCBI Translation Tables |

### Should Tests (Important but not critical)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| S01 | Translate_SequenceShorterThan3_ReturnsEmpty | Incomplete codon | Implementation |
| S02 | FindOrfs_BothStrands_SearchesReverseComplement | Strand option | Wikipedia: ORF |
| S03 | FindOrfs_ForwardOnly_IgnoresReverse | Strand option | Implementation |
| S04 | FindOrfs_OrfResult_HasCorrectPositions | Position tracking | Wikipedia: ORF |
| S05 | TranslateSixFrames_NullInput_ThrowsException | Null handling | Implementation |
| S06 | Translate_VertebrateMito_AGAIsStop | Table 2 specifics | NCBI |
| S07 | Translate_YeastMito_CUUIsTheonine | Table 3 specifics | NCBI |

### Could Tests (Nice to have)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| C01 | Translate_RealSequence_ProducesKnownProtein | Biological validation | ASSUMPTION |
| C02 | FindOrfs_MultipleOrfs_FindsAll | Complex scenario | ASSUMPTION |
| C03 | TranslateSixFrames_EmptySequence_ReturnsEmptyFrames | Edge case | ASSUMPTION |

---

## Existing Test Coverage (Audit)

### Current Test File: TranslatorTests.cs

| Test | Status | Classification |
|------|--------|----------------|
| Translate_SingleCodon_ReturnsSingleAminoAcid | ✓ Exists | M01 - Covered |
| Translate_MultipleCodens_ReturnsProtein | ✓ Exists | M02 - Covered |
| Translate_ToFirstStop_StopsAtStopCodon | ✓ Exists | M03 - Covered |
| Translate_Frame1_ShiftsReading | ✓ Exists | M04 - Covered |
| Translate_Frame2_ShiftsReading | ✓ Exists | M05 - Covered |
| Translate_InvalidFrame_ThrowsException | ✓ Exists | M06 - Covered |
| Translate_EmptySequence_ReturnsEmpty | ✓ Exists | M07 - Covered |
| Translate_NullDna_ThrowsException | ✓ Exists | M08 - Covered |
| Translate_Rna_Works | ✓ Exists | M10 - Covered |
| Translate_RnaToFirstStop_Works | ✓ Exists | M10 - Covered |
| Translate_String_Works | ✓ Exists | M11 - Covered |
| Translate_LowercaseString_Works | ✓ Exists | M12 - Covered |
| Translate_VertebrateMitochondrial_UsesDifferentCode | ✓ Exists | M20 - Covered |
| Translate_YeastMitochondrial_CUU_IsThreonine | ✓ Exists | S07 - Covered |
| TranslateSixFrames_ReturnsAllSixFrames | ✓ Exists | M13 - Covered |
| TranslateSixFrames_Frame1_MatchesDirect | ✓ Exists | M14 - Covered |
| TranslateSixFrames_NegativeFrames_UseReverseComplement | ✓ Exists | M15 - Covered |
| FindOrfs_SimpleOrf_FindsIt | ✓ Exists | M16 - Covered |
| FindOrfs_NoStartCodon_ReturnsEmpty | ✓ Exists | M17 - Covered |
| FindOrfs_ShortOrf_FilteredByMinLength | ✓ Exists | M18 - Covered |
| FindOrfs_RespectMinLength_FindsSmallOrfs | ✓ Exists | M18 - Covered |
| FindOrfs_ForwardOnly_DoesNotSearchReverseStrand | ✓ Exists | S03 - Covered |
| FindOrfs_ResultHasCorrectFrame | ✓ Exists | S04 - Covered |
| FindOrfs_NullDna_ThrowsException | ✓ Exists | M19 - Covered |
| Translate_InsulinBChain_ProducesCorrectProtein | ✓ Exists | C01 - Covered |

---

## Consolidation Plan

### Current State Analysis
The existing test file `TranslatorTests.cs` already has **comprehensive coverage** of the Translator class with 25 tests covering:
- Basic translation (DNA, RNA, string)
- Frame handling
- Stop codon behavior
- Alternative genetic codes
- Six-frame translation
- ORF finding
- Edge cases (empty, null, invalid)

### Gaps Identified
| Missing Test | Priority | Action |
|--------------|----------|--------|
| M09: Translate_NullRna_ThrowsArgumentNullException | Must | Add test |
| TranslateSixFrames_NullInput_ThrowsException | Should | Add test |

### Consolidation Actions
1. **Add**: M09 - Null RNA input handling
2. **Add**: S05 - TranslateSixFrames null handling
3. **Enhance**: Add Test Unit ID header to existing test file
4. **Keep**: All existing tests - well-structured and evidence-aligned

---

## Test Data Sources

### Standard Genetic Code
```
ATG → M (Start/Methionine)
GCT → A (Alanine)
TAA/TAG/TGA → * (Stop)
```

### Alternative Codes
```
Vertebrate Mitochondrial (Table 2):
- AGA/AGG → * (Stop, not Arg)
- UGA → W (Trp, not Stop)

Yeast Mitochondrial (Table 3):  
- CUU/CUC/CUA/CUG → T (Thr, not Leu)
```

---

## Open Questions / Decisions

1. **RESOLVED**: Frame numbering convention - Implementation uses 0, 1, 2 for frame parameter but +1, +2, +3 for TranslateSixFrames keys. This is intentional and matches common conventions.

2. **RESOLVED**: ORF minimum length default - Default is 100 amino acids, which aligns with Wikipedia guidance that ORFs are typically >100 codons for gene prediction.

---

## Validation Criteria

- [ ] All Must tests implemented and passing
- [ ] Zero warnings in test file
- [ ] Test naming follows convention: `Method_Scenario_ExpectedResult`
- [ ] Each test has clear Arrange-Act-Assert structure
- [ ] Evidence sources documented in test file header
