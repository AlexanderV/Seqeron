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
| M13 | VertebrateMito_UGA_IsTryptophan | Table 2 difference | NCBI Table 2 |
| M14 | VertebrateMito_AGA_IsStop | Table 2 difference | NCBI Table 2 |
| M15 | YeastMito_CUU_IsThreonine | Table 3 difference | NCBI Table 3 |
| M16 | BacterialPlastid_AlternativeStartCodons | GUG, UUG as starts | NCBI Table 11 |
| M17 | GetByTableNumber_ValidTable_ReturnsCorrectCode | Factory method | Implementation |
| M18 | GetByTableNumber_InvalidTable_ThrowsException | Error handling | Implementation |

### Should Tests (Quality/Robustness)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| S01 | Standard_Has64Codons | Completeness check |
| S02 | Standard_HasThreeStopCodons | Stop codon count |
| S03 | Standard_HasOneStartCodon | Standard code start |
| S04 | AllCodons_ProduceValidAminoAcids | Output validation |
| S05 | GetCodonsForAminoAcid_Trp_ReturnsOneCodon | Degeneracy: Trp |
| S06 | GetCodonsForAminoAcid_Ser_ReturnsSixCodons | Degeneracy: Ser |
| S07 | VertebrateMito_AUA_IsMethionine | Table 2: AUA=M |
| S08 | VertebrateMito_HasFourStartCodons | Table 2 starts |
| S09 | YeastMito_UGA_IsTryptophan | Table 3: UGA=W |

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
| Standard_HasThreeStopCodons | S | Covered | Stop count |
| Standard_HasOneStartCodon | S | Covered | Start count |
| Translate_AUG_ReturnsMethionine | M | Covered | M03 |
| Translate_DnaCodon_Works | M | Covered | M04 |
| Translate_LowercaseCodon_Works | M | Covered | M05 |
| Translate_StopCodon_ReturnsAsterisk | M | Partial | Only one assertion |
| Translate_InvalidCodonLength_ThrowsException | M | Covered | M12 |
| Translate_AllCodons_ProduceValidAminoAcids | S | Covered | S04 |
| IsStartCodon_AUG_ReturnsTrue | M | Covered | M06 |
| IsStartCodon_ATG_ReturnsTrue | C | Covered | C02 |
| IsStartCodon_Other_ReturnsFalse | M | Covered | M07 |
| IsStopCodon_UAA_ReturnsTrue | M | Partial | Only UAA |
| IsStopCodon_TAA_ReturnsTrue | C | Covered | C03 |
| IsStopCodon_Other_ReturnsFalse | M | Covered | M09 |
| GetCodonsForAminoAcid_Methionine_ReturnsOneCodon | M | Covered | M10 |
| GetCodonsForAminoAcid_Leucine_ReturnsSixCodons | M | Covered | M11 |
| GetCodonsForAminoAcid_Serine_ReturnsSixCodons | S | Covered | S06 |
| GetCodonsForAminoAcid_Tryptophan_ReturnsOneCodon | S | Covered | S05 |
| VertebrateMito_UGA_IsTryptophan | M | Covered | M13 |
| VertebrateMito_AGA_IsStopCodon | M | Covered | M14 |
| VertebrateMito_AUA_IsMethionine | S | Covered | S07 |
| YeastMito_CUU_IsThreonine | M | Covered | M15 |
| BacterialPlastid_HasAlternativeStartCodons | M | Covered | M16 |
| GetByTableNumber_1_ReturnsStandard | M | Covered | M17 |
| GetByTableNumber_2_ReturnsVertebrateMitochondrial | M | Covered | M17 |
| GetByTableNumber_Invalid_ThrowsException | M | Covered | M18 |

### Coverage Gaps

| Gap | Required Test | Priority |
|-----|---------------|----------|
| All 3 stop codons in single assertion | M02 | Must |
| Complete 64-codon coverage | M01 | Must |
| Yeast mito all CUx codons | S | Should |
| Vertebrate mito start codon count | S08 | Should |

---

## Consolidation Plan

1. **Keep existing tests**: GeneticCodeTests.cs has good coverage
2. **Enhance stop codon test**: Test all 3 stop codons comprehensively
3. **Add complete codon coverage test**: Verify all 64 standard codons
4. **Add alternative code completeness tests**: Verify all differences
5. **Remove duplicates**: None identified

---

## Test Data Sources

### Standard Codon Table (NCBI Table 1)
All 64 codon→amino acid mappings from NCBI official reference.

### Alternative Code Differences (NCBI Tables 2, 3, 11)
Key differences documented in Evidence file.

---

## Open Questions / Decisions

1. **Q**: Should we test all 64 codons individually or grouped?
   **A**: Group test for coverage, individual tests for specific behaviors

2. **Q**: Should invalid nucleotides (X, N) be tested?
   **A**: Yes - implementation should throw ArgumentException
