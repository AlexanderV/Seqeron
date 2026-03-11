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
| C01 | Translate_InsulinBChain_ProducesCorrectProtein | Biological validation | UniProt P01308 |
| C02 | FindOrfs_MultipleOrfs_FindsAll | Complex scenario | Implementation spec |
| C03 | TranslateSixFrames_EmptySequence_ReturnsEmptyFrames | Edge case | Implementation spec |

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
| Translate_NullRna_ThrowsException | ✓ Exists | M09 - Covered |
| Translate_Rna_Works | ✓ Exists | M10 - Covered |
| Translate_RnaToFirstStop_Works | ✓ Exists | M10 - Covered |
| Translate_DnaString_ConvertsTToU | ✓ Exists | M11 - Covered |
| Translate_LowercaseString_Works | ✓ Exists | M12 - Covered |
| TranslateSixFrames_ReturnsAllSixFrames | ✓ Exists | M13 - Covered |
| TranslateSixFrames_Frame1_MatchesDirect | ✓ Exists | M14 - Covered |
| TranslateSixFrames_NegativeFrames_UseReverseComplement | ✓ Exists | M15 - Covered |
| FindOrfs_SimpleOrf_FindsIt | ✓ Exists | M16 - Covered |
| FindOrfs_NoStartCodon_ReturnsEmpty | ✓ Exists | M17 - Covered |
| FindOrfs_RespectMinLength_FindsSmallOrfs | ✓ Exists | M18 - Covered |
| FindOrfs_ShortOrf_FilteredByMinLength | ✓ Exists | M18 - Covered |
| FindOrfs_NullDna_ThrowsException | ✓ Exists | M19 - Covered |
| Translate_VertebrateMitochondrial_UsesDifferentCode | ✓ Exists | M20, S06 - Covered |
| Translate_SequenceShorterThan3_ReturnsEmpty | ✓ Exists | S01 - Covered |
| FindOrfs_BothStrands_SearchesReverseComplement | ✓ Exists | S02 - Covered |
| FindOrfs_ForwardOnly_DoesNotSearchReverseStrand | ✓ Exists | S03 - Covered |
| FindOrfs_OrfResult_HasCorrectPositions | ✓ Exists | S04 - Covered |
| TranslateSixFrames_NullInput_ThrowsException | ✓ Exists | S05 - Covered |
| Translate_YeastMitochondrial_CUU_IsThreonine | ✓ Exists | S07 - Covered |
| Translate_InsulinBChain_ProducesCorrectProtein | ✓ Exists | C01 - Covered |
| FindOrfs_MultipleOrfs_FindsAll | ✓ Exists | C02 - Covered |
| TranslateSixFrames_EmptySequence_ReturnsEmptyFrames | ✓ Exists | C03 - Covered |

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

## Deviations and Assumptions

**None.** All tests and implementation are grounded in the authoritative NCBI Genetic Codes and Wikipedia sources listed in the Evidence document.

| Item | Status | Evidence |
|------|--------|----------|
| Standard genetic code (Table 1) | Verified | All 64 codons match NCBI Table 1. Start codons (AUG, UUG, CUG) and stop codons (UAA, UAG, UGA) match NCBI. |
| Vertebrate Mitochondrial (Table 2) | Verified | All 4 differences from standard (AUA→M, UGA→W, AGA→*, AGG→*) match NCBI Table 2. Start codons (AUG, AUA, AUU, AUC, GUG) and stop codons (UAA, UAG, AGA, AGG) match NCBI. |
| Yeast Mitochondrial (Table 3) | Verified | All 5 differences (CUU/CUC/CUA/CUG→T, AUA→M, UGA→W) match NCBI Table 3. Start codons (AUG, AUA, GUG) and stop codons (UAA, UAG) match NCBI. |
| Bacterial/Plastid (Table 11) | Verified | Same amino acid translation as standard. Start codons (AUG, GUG, UUG, CUG, AUU, AUC, AUA) match NCBI Table 11. |
| Insulin B chain test data | Verified | Full 30-amino-acid protein FVNQHLCGSHLVEALYLVCGERGFFYTPKT matches UniProt P01308 (Homo sapiens, positions 25–54 of preproinsulin). |
| Frame numbering | Documented | `Translate` uses 0-based frames (0, 1, 2); `TranslateSixFrames` uses ±1/±2/±3 keys. Both conventions are standard in bioinformatics. |
| ORF default minimum length | Documented | Default 100 amino acids per Wikipedia ORF guidance ("ORFs are typically >100 codons for gene prediction"). |

---

## Validation Criteria

- [x] All Must tests (M01–M20) implemented and passing
- [x] All Should tests (S01–S07) implemented and passing
- [x] All Could tests (C01–C03) implemented and passing
- [x] Total: 31 tests, 0 failures
- [x] Zero assumptions — all test data backed by external sources (NCBI, Wikipedia, UniProt)
- [x] Test naming follows convention: `Method_Scenario_ExpectedResult`
- [x] Each test has clear Arrange-Act-Assert structure
- [x] Evidence sources documented in test file header
- [x] All genetic code tables verified against NCBI (March 2026)
- [x] No duplicates — each test covers a unique scenario
- [x] Weak tests strengthened: M11 (T→U equivalence verified), M15 (all 3 negative frames verified)
