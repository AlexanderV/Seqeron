# TestSpec: PRIMER-DESIGN-001

## Primer Pair Design

**Test Unit ID:** PRIMER-DESIGN-001
**Area:** MolTools
**Status:** ☑ Complete
**Last Updated:** 2026-01-22

---

## Evidence Summary

### Authoritative Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Primer (molecular biology) | Encyclopedia | 18-24 bp length, 40-60% GC, Tm 50-60°C, primer pairs within 5°C |
| Addgene: How to Design a Primer | Protocol Guide | Length 18-24, GC 40-60%, Tm 50-60°C, avoid complementary regions |
| Primer3 Manual | Software Documentation | Min/Max/Optimal sizes, GC content, Tm formulas, hairpin/dimer detection |
| SantaLucia (1998) PNAS 95:1460-65 | Research Paper | Nearest-neighbor thermodynamics for Tm calculation |

### Standard Primer Design Parameters

From authoritative sources (Primer3, Addgene):

| Parameter | Default/Recommended | Source |
|-----------|-------------------|--------|
| Length (Min) | 18 bp | Primer3, Wikipedia |
| Length (Max) | 25-27 bp | Primer3 |
| Length (Optimal) | 20 bp | Primer3 |
| GC Content (Min) | 20-40% | Primer3: 20%, Addgene: 40% |
| GC Content (Max) | 60-80% | Addgene: 60%, Primer3: 80% |
| Tm (Min) | 50-57°C | Addgene: 50°C, Primer3: 57°C |
| Tm (Max) | 60-65°C | Addgene: 60°C, Primer3: 63°C |
| Pair Tm Difference | ≤ 5°C | Wikipedia, Addgene |
| Homopolymer Max | 4-5 | Primer3 (PRIMER_MAX_POLY_X = 5) |

### Key Design Principles

1. **3' End Stability**: GC clamp at 3' end beneficial but not excessive (max 5 GC in last 5 bases)
2. **Hairpin Avoidance**: Primers should not form stable secondary structures
3. **Primer-Dimer Prevention**: 3' ends should not be complementary between primer pairs
4. **Product Size**: Forward and reverse primers should amplify target region

---

## Methods Under Test

| Method | Class | Type | Complexity |
|--------|-------|------|------------|
| `DesignPrimers(template, start, end, params)` | PrimerDesigner | Canonical | O(n²) |
| `EvaluatePrimer(seq, pos, isForward, params)` | PrimerDesigner | Helper | O(m²) |
| `GeneratePrimerCandidates(template, region)` | PrimerDesigner | Helper | O(n×m) |

---

## Test Requirements

### MUST Tests (Required)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | DesignPrimers returns valid primer pair for suitable template | Core functionality | Primer3, Addgene |
| M2 | DesignPrimers throws on invalid target range (start >= end) | Input validation | Implementation |
| M3 | DesignPrimers throws on out-of-bounds target | Input validation | Implementation |
| M4 | Forward primer position is upstream of target | Algorithm invariant | Primer3 |
| M5 | Reverse primer position is downstream of target | Algorithm invariant | Primer3 |
| M6 | Primer pair Tm difference ≤ 5°C when valid | Standard requirement | Wikipedia, Addgene |
| M7 | EvaluatePrimer validates length constraints | Primer3 defaults | Primer3: 18-27 |
| M8 | EvaluatePrimer validates GC content constraints | Standard | Primer3: 20-80%, Impl: 40-60% |
| M9 | EvaluatePrimer validates Tm constraints | Standard | Primer3: 57-63°C, Impl: 55-65°C |
| M10 | EvaluatePrimer detects homopolymer runs | Avoid repeats | Primer3: max 5 |
| M11 | EvaluatePrimer detects dinucleotide repeats | Avoid repeats | Primer3, Wikipedia |
| M12 | GeneratePrimerCandidates generates forward primers correctly | Correctness | Implementation |
| M13 | GeneratePrimerCandidates generates reverse primers as reverse complement | Correctness | Implementation |

### SHOULD Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| S1 | Primer pair avoids primer-dimer formation | Best practice |
| S2 | Product size is correctly calculated | Usability |
| S3 | Custom parameters are respected | Flexibility |
| S4 | Score calculation rewards optimal properties | Quality ranking |
| S5 | Return failure message when no valid primers found | Error handling |

### COULD Tests (Nice to Have)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| C1 | Performance on long templates | Efficiency |
| C2 | Multiple overlapping candidates generated | Completeness |

---

## Audit of Existing Tests

### File: PrimerDesignerTests.cs

| Test Method | Coverage | Status |
|-------------|----------|--------|
| `DesignPrimers_ValidTarget_ReturnsResult` | M1 (partial) | Weak - doesn't assert primer validity |
| `DesignPrimers_InvalidTarget_ThrowsException` | M2 | Covered |
| `DesignPrimers_TargetOutOfRange_ThrowsException` | M3 | Covered |
| `GeneratePrimerCandidates_ReturnsMultipleCandidates` | M12 (partial) | Weak |
| `GeneratePrimerCandidates_Forward_HasCorrectOrientation` | M12 | Covered |
| `GeneratePrimerCandidates_Reverse_HasCorrectOrientation` | M13 | Covered |
| `EvaluatePrimer_GoodPrimer_IsValid` | General | Covered |
| `EvaluatePrimer_TooShort_NotValid` | M7 | Covered |
| `EvaluatePrimer_TooLongHomopolymer_NotValid` | M10 | Covered |
| `EvaluatePrimer_CalculatesScore` | S4 (partial) | Weak |
| `EvaluatePrimer_CustomParameters_AppliesCorrectly` | S3 | Covered |
| `EvaluatePrimer_TypicalGoodPrimer_PassesAllChecks` | General | Covered |
| `EvaluatePrimer_ProblematicPrimer_DetectsIssues` | General | Covered |

### Missing Coverage

- M4: Forward primer upstream of target assertion
- M5: Reverse primer downstream of target assertion
- M6: Primer pair Tm difference validation
- M8: GC content constraint validation
- M9: Tm constraint validation
- M11: Dinucleotide repeat validation (partially tested)
- S1: Primer-dimer avoidance in pair design
- S2: Product size calculation
- S5: Failure message on no valid primers

---

## Consolidation Plan

### Canonical Test File

`PrimerDesigner_PrimerDesign_Tests.cs` - New file focused on PRIMER-DESIGN-001

### Tests to Add

1. Forward/Reverse primer position validation
2. Tm difference validation for valid pairs
3. GC content constraint tests
4. Tm constraint tests
5. Product size verification
6. No-valid-primers failure case
7. Comprehensive end-to-end design test

### Tests to Move from PrimerDesignerTests.cs

Tests covering DesignPrimers, EvaluatePrimer (for design), and GeneratePrimerCandidates should be consolidated into the new canonical file, leaving only GC content, homopolymer, and dinucleotide helper tests in PrimerDesignerTests.cs.

---

## Open Questions

None - behavior is well-defined by implementation and aligns with authoritative sources.

---

## Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Implementation's 40-60% GC range is stricter than Primer3's 20-80% | More conservative for PCR success |
| A2 | Implementation's 55-65°C Tm range aligns with practical PCR conditions | Common protocol range |
