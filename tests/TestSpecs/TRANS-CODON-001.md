# TestSpec: TRANS-CODON-001 - Codon Translation

**Test Unit ID:** TRANS-CODON-001  
**Area:** Translation  
**Canonical Class:** `GeneticCode`  
**Date:** 2026-02-04

---

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `Translate(codon)` | GeneticCode | Canonical | Deep |
| `IsStartCodon(codon)` | GeneticCode | Canonical | Deep |
| `IsStopCodon(codon)` | GeneticCode | Canonical | Deep |
| `GetCodonsForAminoAcid(aa)` | GeneticCode | Canonical | Deep |
| `GetByTableNumber(n)` | GeneticCode | Factory | Moderate |

---

## Test Classification

### Must Tests (Evidence-backed)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M01 | Translate_AllStandardCodons_CorrectAminoAcid | Complete coverage of 64 codons | NCBI Table 1 |
| M02 | Translate_StopCodons_ReturnsAsterisk | UAA, UAG, UGA → '*' | Wikipedia: Stop codon |
| M03 | Translate_StartCodon_ReturnsMethionine | AUG → 'M' | Wikipedia: Start codon |
| M04 | Translate_DnaCodon_ConvertedToRna | ATG ≡ AUG | Implementation spec |
| M05 | Translate_LowercaseCodon_CaseInsensitive | aug = AUG | Implementation spec |
| M06 | IsStartCodon_AUG_ReturnsTrue | Primary start codon | NCBI, Wikipedia |
| M07 | IsStartCodon_NonStart_ReturnsFalse | UUU is not start | NCBI Table 1 |
| M08 | IsStopCodon_AllThreeStops_ReturnsTrue | UAA, UAG, UGA | Wikipedia: Stop codon |
| M09 | IsStopCodon_NonStop_ReturnsFalse | AUG is not stop | NCBI Table 1 |
| M10 | GetCodonsForAminoAcid_Met_ReturnsOneCodon | M has only AUG | NCBI degeneracy |
| M11 | GetCodonsForAminoAcid_Leu_ReturnsSixCodons | L has 6 codons | NCBI degeneracy |
| M12 | Translate_InvalidLength_ThrowsException | Codon must be 3 chars | Definition |
| M19 | Translate_InvalidNucleotide_ThrowsException | Non-ACGTU nucleotides rejected | NCBI definition |
| M13 | VertebrateMito_UGA_IsTryptophan | Table 2 difference | NCBI Table 2 |
| M14 | VertebrateMito_AGA_IsStop | Table 2 difference | NCBI Table 2 |
| M15 | YeastMito_CUU_IsThreonine | Table 3 difference | NCBI Table 3 |
| M16 | BacterialPlastid_HasSevenStartCodons | Table 11: 7 starts per NCBI Starts line | NCBI Table 11 |
| M17 | GetByTableNumber_ValidTable_ReturnsCorrectCode | Factory method | Implementation |
| M18 | GetByTableNumber_InvalidTable_ThrowsException | Error handling | Implementation |

### Should Tests (Quality/Robustness)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| S01 | Standard_Has64Codons | Completeness check |
| S02 | Standard_HasThreeStopCodons | Stop codon count |
| S03 | Standard_HasThreeStartCodons | Standard code: 3 starts (AUG, UUG, CUG) | NCBI Table 1 Starts line |
| S04 | AllCodons_ProduceValidAminoAcids | Output validation |
| S05 | GetCodonsForAminoAcid_Trp_ReturnsOneCodon | Degeneracy: Trp |
| S06 | GetCodonsForAminoAcid_Ser_ReturnsSixCodons | Degeneracy: Ser |
| S07 | VertebrateMito_AUA_IsMethionine | Table 2: AUA=M |
| S08 | VertebrateMito_HasFiveStartCodons | Table 2: 5 starts (AUG,AUA,AUU,AUC,GUG) |
| S09 | YeastMito_UGA_IsTryptophan | Table 3: UGA=W |
| S10 | YeastMito_AUA_IsMethionine | Table 3: AUA=M | NCBI Table 3 |

### Could Tests (Enhancement)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C01 | Translate_MixedCaseCodon_Works | Robustness |
| C02 | IsStartCodon_DnaFormat_Works | ATG as start |
| C03 | IsStopCodon_DnaFormat_Works | TAA/TAG/TGA |

---

## Existing Test Coverage (Audit)

### GeneticCodeTests.cs Analysis

| Test | Category | Status | Notes |
|------|----------|--------|-------|
| Standard_HasCorrectName | S | Covered | Metadata |
| Standard_Has64Codons | S | Covered | Completeness |
| Standard_HasThreeStopCodons | S | Covered | Exact set: UAA, UAG, UGA |
| Standard_HasThreeStartCodons | S | Covered | Exact set: AUG, UUG, CUG |
| Translate_AUG_ReturnsMethionine | M | Covered | M03 |
| Translate_DnaCodon_Works | M | Covered | M04 |
| Translate_LowercaseCodon_Works | M | Covered | M05 |
| Translate_AllStopCodons_ReturnsAsterisk | M | Covered | M02: UAA, UAG, UGA |
| Translate_InvalidCodonLength_ThrowsException | M | Covered | M12 |
| Translate_EmptyCodon_ThrowsException | M | Covered | M12 |
| Translate_NullCodon_ThrowsException | M | Covered | M12 |
| Translate_AllCodons_ProduceValidAminoAcids | S | Covered | S04 |
| Translate_CompleteStandardCodonTable_MatchesNcbi | M | Covered | M01: all 64 codons |
| Translate_InvalidNucleotide_ThrowsException | M | Covered | M19: XYZ, ANN, 12G |
| Translate_MixedCaseCodon_Works | C | Covered | C01 |
| IsStartCodon_AUG_ReturnsTrue | M | Covered | M06 |
| IsStartCodon_ATG_ReturnsTrue | C | Covered | C02 |
| IsStartCodon_NonStartCodon_ReturnsFalse | M | Covered | M07 |
| IsStartCodon_InvalidInput_ReturnsFalse | M | Covered | Defensive: short/empty/null |
| IsStopCodon_AllThreeStandard_ReturnsTrue | M | Covered | M08: UAA, UAG, UGA |
| IsStopCodon_DnaFormat_ReturnsTrue | C | Covered | C03: TAA, TAG, TGA |
| IsStopCodon_NonStopCodon_ReturnsFalse | M | Covered | M09 |
| GetCodonsForAminoAcid_Methionine_ReturnsOneCodon | M | Covered | M10: exact {AUG} |
| GetCodonsForAminoAcid_Leucine_ReturnsSixCodons | M | Covered | M11: exact 6-codon set |
| GetCodonsForAminoAcid_Serine_ReturnsSixCodons | S | Covered | S06: exact {UCU,UCC,UCA,UCG,AGU,AGC} |
| GetCodonsForAminoAcid_Arginine_ReturnsSixCodons | S | Covered | Exact {CGU,CGC,CGA,CGG,AGA,AGG} |
| GetCodonsForAminoAcid_Tryptophan_ReturnsOneCodon | S | Covered | S05: exact {UGG} |
| GetCodonsForAminoAcid_Isoleucine_ReturnsThreeCodons | S | Covered | Exact {AUU,AUC,AUA} |
| GetCodonsForAminoAcid_Stop_ReturnsThreeCodons | S | Covered | Exact {UAA,UAG,UGA} |
| GetCodonsForAminoAcid_LowercaseInput_Works | C | Covered | Case-insensitive |
| VertebrateMitochondrial_UGA_IsTryptophan | M | Covered | M13 |
| VertebrateMitochondrial_AGA_IsStopCodon | M | Covered | M14 |
| VertebrateMitochondrial_AGG_IsStopCodon | M | Covered | Table 2 |
| VertebrateMitochondrial_AUA_IsMethionine | S | Covered | S07 |
| VertebrateMitochondrial_HasFiveStartCodons | S | Covered | S08: exact 5-codon set |
| VertebrateMitochondrial_HasFourStopCodons | S | Covered | Table 2: UAA,UAG,AGA,AGG |
| YeastMitochondrial_AllCUxCodons_AreThreonine | M | Covered | M15: CUU,CUC,CUA,CUG=T |
| YeastMitochondrial_AUA_IsMethionine | S | Covered | S10: Table 3 AUA=M |
| YeastMitochondrial_UGA_IsTryptophan | S | Covered | S09 |
| YeastMitochondrial_HasTwoStopCodons | S | Covered | Table 3: UAA,UAG |
| YeastMitochondrial_HasThreeStartCodons | S | Covered | Table 3: AUG,AUA,GUG |
| BacterialPlastid_HasSevenStartCodons | M | Covered | M16: exact 7-codon set |
| BacterialPlastid_CodonTable_SameAsStandard | S | Covered | All 64 codons compared |
| GetByTableNumber_1_ReturnsStandard | M | Covered | M17 |
| GetByTableNumber_2_ReturnsVertebrateMitochondrial | M | Covered | M17 |
| GetByTableNumber_3_ReturnsYeastMitochondrial | M | Covered | M17 |
| GetByTableNumber_11_ReturnsBacterialPlastid | M | Covered | M17 |
| GetByTableNumber_AllSupported_ReturnsValidCodes | S | Covered | M17 |
| GetByTableNumber_Invalid_ThrowsException | M | Covered | M18 |

### Coverage Gaps (All Closed)

All coverage gaps resolved:
- ❌ Missing: `Translate_InvalidNucleotide_ThrowsException` (M19) added
- ❌ Missing: `YeastMitochondrial_AUA_IsMethionine` (S10) added
- ⚠ Weak: `GetCodonsForAminoAcid_Serine` strengthened with exact 6-codon set
- ⚠ Weak: `GetCodonsForAminoAcid_Isoleucine` strengthened with exact 3-codon set
- ⚠ Weak: `GetCodonsForAminoAcid_Arginine` strengthened with all 6 codons
- ⚠ Weak: `BacterialPlastid_CodonTable_SameAsStandard` now compares all 64 codons
- 🔁 Duplicate: `IsStopCodon_UAA_ReturnsTrue` removed (covered by `AllThreeStandard`)
- 🔁 Duplicate: `YeastMitochondrial_CUU_IsThreonine` removed (covered by `AllCUxCodons`)

---

## Consolidation Plan

All tests consolidated. No duplicates. Full coverage achieved (47 tests).

---

## Test Data Sources

### Standard Codon Table (NCBI Table 1)
All 64 codon→amino acid mappings from NCBI official reference.

### Alternative Code Differences (NCBI Tables 2, 3, 11)
Key differences documented in Evidence file.

---

## Deviations and Assumptions

**None.** Implementation and tests match NCBI translation tables exactly.

All start/stop codon sets are derived directly from NCBI `Starts` and `AAs` strings
for Tables 1, 2, 3, and 11 (NCBI last updated Sep. 23, 2024).

---

## Open Questions / Decisions

1. **Q**: Should we test all 64 codons individually or grouped?
   **A**: Group test for coverage, individual tests for specific behaviors

2. **Q**: Should invalid nucleotides (X, N) be tested?
   **A**: Yes — `Translate_InvalidNucleotide_ThrowsException` (M19) covers this
