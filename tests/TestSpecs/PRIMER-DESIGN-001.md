# TestSpec: PRIMER-DESIGN-001

## Primer Pair Design

**Test Unit ID:** PRIMER-DESIGN-001
**Area:** MolTools
**Status:** ☑ Complete
**Last Updated:** 2026-03-04
**Total Tests:** 149 (canonical + smoke + mutation-killing)

---

## Evidence Summary

### Authoritative Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia: Primer (molecular biology)](https://en.wikipedia.org/wiki/Primer_(molecular_biology)) | Encyclopedia | 18-24 bp length, 40-60% GC, Tm 50-60°C, primer pairs within 5°C |
| [Addgene: How to Design a Primer](https://www.addgene.org/protocols/primer-design/) | Protocol Guide | Length 18-24, GC 40-60%, Tm 50-60°C, pairs within 5°C, avoid complementary regions |
| [Primer3 Manual (v2.6.1)](https://primer3.org/manual.html) | Software Documentation | PRIMER_MIN_SIZE=18, PRIMER_MAX_SIZE=27, PRIMER_OPT_SIZE=20, PRIMER_MIN_TM=57, PRIMER_OPT_TM=60, PRIMER_MAX_TM=63, PRIMER_MIN_GC=20, PRIMER_MAX_GC=80, PRIMER_MAX_POLY_X=5, PRIMER_PAIR_MAX_DIFF_TM=100.0 |
| SantaLucia (1998) PNAS 95:1460-65 | Research Paper | Nearest-neighbor thermodynamics for Tm calculation |

### Implementation Parameters vs Sources

| Parameter | Implementation | Authoritative Source | Justification |
|-----------|---------------|---------------------|---------------|
| Length (Min) | 18 bp | Primer3: 18 | Exact match |
| Length (Max) | 25 bp | Primer3: 27, Addgene: 24 | Practical middle ground; within both ranges |
| Length (Optimal) | 20 bp | Primer3: 20 | Exact match |
| GC Content (Min) | 40% | Addgene: 40% | Exact match (Addgene) |
| GC Content (Max) | 60% | Addgene: 60% | Exact match (Addgene) |
| Tm (Min) | 57°C | Primer3: 57°C | Exact match |
| Tm (Max) | 63°C | Primer3: 63°C | Exact match |
| Tm (Optimal) | 60°C | Primer3: 60°C | Exact match |
| Pair Tm Difference | ≤ 5°C | Wikipedia, Addgene | Exact match (Primer3 default PRIMER_PAIR_MAX_DIFF_TM=100.0 is unlimited; 5°C is the standard lab guideline) |
| Homopolymer Max | 4 | Primer3: 5 | Stricter than Primer3; conservative choice |

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
| M7 | EvaluatePrimer validates length constraints (18-25 bp) | Primer3 defaults | Primer3: 18-27 |
| M8 | EvaluatePrimer validates GC content constraints (40-60%) | Addgene standard | Addgene: 40-60% |
| M9 | EvaluatePrimer validates Tm constraints (57-63°C) | Primer3 defaults | Primer3: 57-63°C |
| M10 | EvaluatePrimer detects homopolymer runs (max 4) | Avoid repeats | Primer3: max 5 |
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

## Canonical Test File

`PrimerDesigner_PrimerDesign_Tests.cs`

All PRIMER-DESIGN-001 tests consolidated. Smoke tests and mutation-killing tests remain in `PrimerDesignerTests.cs`.

---

## Coverage Classification

Applied systematic coverage classification (2026-03-04):

### MUST Tests

| ID | Status | Test Method(s) | Notes |
|----|--------|----------------|-------|
| M1 | ✅ | `DesignPrimers_ValidTemplate_ForwardIsUpstreamOfTarget`, `_ForwardWithinSearchRegion` | Asserts `result.IsValid == true` — non-vacuous |
| M2 | ✅ | `DesignPrimers_TargetEndBeforeStart_ThrowsArgumentException` | Exact exception type |
| M3 | ✅ | `DesignPrimers_TargetBeyondTemplate_ThrowsArgumentException`, `_NegativeCoordinates_ThrowsArgumentException` | Two edge cases |
| M4 | ✅ | `DesignPrimers_ValidTemplate_ForwardIsUpstreamOfTarget` | Forward.Position < targetStart |
| M5 | ✅ | `DesignPrimers_ValidTemplate_ReverseIsDownstreamOfTarget` | Reverse.Position >= targetEnd |
| M6 | ✅ | `DesignPrimers_PrimerPair_TmDifferenceWithin5Degrees` | Exact ≤5°C assertion |
| M7 | ✅ | `DesignPrimers_Primers_HaveLengthWithinRange`, `EvaluatePrimer_LengthOutsideRange_ReportsIssue(17,26)` | Positive + boundary |
| M8 | ✅ | `DesignPrimers_Primers_HaveGcContentWithinRange`, `EvaluatePrimer_GcOutsideRange_ReportsIssue(100%,0%)` | Positive + boundary |
| M9 | ✅ | `DesignPrimers_Primers_HaveTmWithinRange` + 4 mutation-killing boundary tests in `PrimerDesignerTests.cs` | Exact range + boundary mutations killed |
| M10 | ✅ | `DesignPrimers_Primers_NoExcessiveHomopolymers`, `EvaluatePrimer_ExcessiveHomopolymer_ReportsIssue` | Positive + negative |
| M11 | ✅ | `EvaluatePrimer_ExcessiveDinucleotideRepeats_ReportsIssue` | Exact repeat count + issue detection |
| M12 | ✅ | `GeneratePrimerCandidates_ReturnsMultipleValidCandidates`, `_AllCandidatesHaveValidLength` | Count + length range |
| M13 | ✅ | `GeneratePrimerCandidates_Reverse_SequenceIsReverseComplement` | Verifies actual revcomp of template substring |

### SHOULD Tests

| ID | Status | Test Method(s) | Notes |
|----|--------|----------------|-------|
| S1 | ✅ | `DesignPrimers_PrimerPair_NoPrimerDimerFormation`, `HasPrimerDimer_ComplementaryPrimers_ReturnsTrue` | Positive (no dimer in valid pair) + negative (engineered dimer) |
| S2 | ✅ | `DesignPrimers_ValidResult_ProductSizeCorrect` | Exact formula: Reverse.Position + Reverse.Length - Forward.Position |
| S3 | ✅ | `DesignPrimers_CustomParameters_AppliesLengthRange`, `GeneratePrimerCandidates_CustomParameters_AppliesLengthRange` | Custom length range (22-28) + (20-22) |
| S4 | ✅ | `EvaluatePrimer_OptimalPrimer_HasHighScore`, `_SuboptimalLength_ScoreVaries` | Score > 50 for optimal; optimal >= shorter |
| S5 | ✅ | `DesignPrimers_HomopolymerRichTemplate_MayReturnInvalid`, `_VeryShortTemplate_ThrowsArgumentException` | Failure message + exception |

### COULD Tests

| ID | Status | Test Method(s) | Notes |
|----|--------|----------------|-------|
| C1 | ✅ | `DesignPrimers_LongTemplate_CompletesWithinTimeout` | 10kb template, < 5s timeout |
| C2 | ✅ | `GeneratePrimerCandidates_LargeRegion_ReturnsMultipleCandidates` | Count > 1 from 80bp region |

### Additional Tests

| Category | Test Method(s) | File |
|----------|----------------|------|
| Hairpin detection | `EvaluatePrimer_SelfComplementary_DetectsHairpin`, `HasHairpinPotential_LongSelfComplementary_ReturnsTrue` | Canonical |
| 3' stability | `Calculate3PrimeStability_GCRich3Prime_MoreNegative` | Canonical |
| Edge cases | `DesignPrimers_NullTemplate_ThrowsException`, `EvaluatePrimer_EmptySequence_HandledGracefully`, `GeneratePrimerCandidates_EmptyRegion_ReturnsEmpty` | Canonical |
| Default params snapshot | `DefaultParameters_HasReasonableValues` | Smoke |
| GC content helpers | `CalculateGcContent_AllGC_Returns100`, `_NoGC_Returns0`, `_HalfGC_Returns50`, `_EmptySequence_Returns0` | Smoke |
| Mutation-killing (Tm) | `EvaluatePrimer_TmOnlyBelowMin_FlagsTmIssue`, `_TmOnlyAboveMax_FlagsTmIssue`, `_TmExactlyAtMinTm_NoTmIssue`, `_TmExactlyAtMaxTm_NoTmIssue` | Smoke |
| Mutation-killing (Hairpin) | `HasHairpinPotential_NullSequence_ReturnsFalse`, `_EmptySequence_ReturnsFalse`, `_SequenceExactlyAtThreshold_DoesNotReturnEarly` | Smoke |
| Cross-reference smoke | `CalculateMeltingTemperature_SmokeTest_ReturnsValidValue`, `FindLongestHomopolymer_SmokeTest_ReturnsValidValue`, `FindLongestDinucleotideRepeat_SmokeTest_ReturnsValidValue`, `HasHairpinPotential_SmokeTest_ReturnsExpectedValue`, `HasPrimerDimer_SmokeTest_ReturnsExpectedValue`, `Calculate3PrimeStability_SmokeTest_ReturnsNegativeValue` | Smoke |

### Summary

- **Missing:** 0
- **Weak:** 0 (all vacuous `if (result.IsValid)` guards removed; all `Is.True.Or.False` assertions replaced)
- **Duplicate:** 0 (13 duplicates removed from `PrimerDesignerTests.cs`)
- **Covered:** 20/20 (M1-M13 + S1-S5 + C1-C2)

---

## Deviations and Assumptions

**Assumptions:** None.

**Intentional deviations from defaults (documented in Parameters table above):**

- **Max Length**: 25 bp — between Primer3 (27) and Addgene (24); practical middle ground.
- **Homopolymer Max**: 4 — stricter than Primer3 default (5); conservative choice.
- **GC Content**: 40-60% — follows Addgene guideline; stricter than Primer3 (20-80%).
- **Pair Tm Difference**: ≤ 5°C — follows Addgene/Wikipedia; Primer3 default (100.0°C) is unlimited.
