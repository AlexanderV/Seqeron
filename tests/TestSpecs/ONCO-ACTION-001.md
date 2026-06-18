# Test Specification: ONCO-ACTION-001

**Test Unit ID:** ONCO-ACTION-001
**Area:** Oncology
**Algorithm:** Clinical Actionability Assessment (OncoKB Therapeutic Levels of Evidence)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Chakravarty D et al. (2017). OncoKB: A Precision Oncology Knowledge Base. JCO PO. | 1 | https://doi.org/10.1200/PO.17.00011 | 2026-06-15 |
| 2 | OncoKB Therapeutic Levels of Evidence (V2) PDF | 3 | https://www.oncokb.org/content/files/levelOfEvidence/V2/LevelsOfEvidence.pdf | 2026-06-15 |
| 3 | OncoKB Curation SOP v3 PDF | 3 | https://sop.oncokb.org/static/sop/OncoKB_Curation_Standard_Operating_Procedure_v3.pdf | 2026-06-15 |
| 4 | oncokb-annotator README (HIGHEST_LEVEL columns) | 3 | https://github.com/oncokb/oncokb-annotator | 2026-06-15 |

### 1.2 Key Evidence Points

1. Seven therapeutic levels are defined: 1, 2, 3A, 3B, 4 (sensitivity) and R1, R2 (resistance) — OncoKB Levels PDF (Source 2).
2. Combined actionability order (highest→lowest): **R1 > 1 > 2 > 3A > 3B > 4 > R2** — oncokb-annotator README HIGHEST_LEVEL (Source 4).
3. Sensitive-only order: **1 > 2 > 3A > 3B > 4** — README HIGHEST_SENSITIVE_LEVEL (Source 4).
4. Resistance-only order: **R1 > R2** — README HIGHEST_RESISTANCE_LEVEL (Source 4).
5. Levels 1, 2 (and R1) are "standard care" / highest; 3A, 3B, 4 (and R2) are investigational/hypothetical — SOP v3 (Source 3).
6. The highest level for a variant is the maximum, under the order above, over all its leveled drug associations — README (Source 4).

### 1.3 Documented Corner Cases

- Variant with no leveled drug association → no highest level (annotator leaves HIGHEST_LEVEL empty); modeled here as `NotActionable` (Source 4 + ASSUMPTION A1).
- A variant may carry both a sensitive and a resistance association; the two axes are reported separately (SOP v3, README).

### 1.4 Known Failure Modes / Pitfalls

1. Mis-ordering 3A vs 3B — the system was refined so 3A (in-indication clinical evidence) outranks 3B (other-indication standard care) — SOP v3 (Source 3).
2. Treating R1 as low because it is "R" — R1 is above Level 1 in the combined order — README (Source 4).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `AssessActionability(IEnumerable<VariantActionabilityInput>)` | OncologyAnalyzer | Canonical | Per-variant highest combined/sensitive/resistance level + NotActionable. |
| `ClassifyActionabilityLevel(VariantActionabilityInput)` | OncologyAnalyzer | Canonical | Returns the single highest combined OncoKB level (or NotActionable) for one variant. |
| `GetTherapyRecommendations(VariantActionabilityInput)` | OncologyAnalyzer | Delegate | Returns caller-supplied therapy associations ordered by descending level. |
| `CompareLevels(OncoKbLevel, OncoKbLevel)` | OncologyAnalyzer | Internal | Ordering comparator; tested indirectly + smoke. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | The highest level returned for a variant is the maximum of its associations' levels under R1 > 1 > 2 > 3A > 3B > 4 > R2. | Yes | Source 4 (README HIGHEST_LEVEL) |
| INV-2 | Highest sensitive level uses 1 > 2 > 3A > 3B > 4 and ignores R1/R2. | Yes | Source 4 (HIGHEST_SENSITIVE_LEVEL) |
| INV-3 | Highest resistance level uses R1 > R2 and ignores sensitivity levels. | Yes | Source 4 (HIGHEST_RESISTANCE_LEVEL) |
| INV-4 | A variant with zero associations yields NotActionable (no level on any axis). | Yes | Source 4 + ASSUMPTION A1 |
| INV-5 | Output has exactly one assessment per input variant, in input order. | Yes | Library convention (mirrors `AnnotateCancerVariants`) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Combined order full chain | Variant whose associations are each single level, compared pairwise | R1 > 1 > 2 > 3A > 3B > 4 > R2 holds for all adjacent pairs | Source 4 HIGHEST_LEVEL |
| M2 | Highest combined of {1, R1} | Both standard-care sensitivity and standard-care resistance | Combined = R1 (R1 outranks 1) | Source 4 |
| M3 | Highest combined of {4, R2} | Hypothetical sensitivity + investigational resistance | Combined = Level4 (4 outranks R2) | Source 4 |
| M4 | Highest sensitive of {2, 3A} | Two sensitivity associations | Sensitive = Level2 | Source 4 HIGHEST_SENSITIVE_LEVEL |
| M5 | Highest sensitive of {3A, 3B, 4} | Three investigational sensitivity associations | Sensitive = Level3A | Source 4 |
| M6 | Highest resistance of {R1, R2} | Two resistance associations | Resistance = LevelR1 | Source 4 HIGHEST_RESISTANCE_LEVEL |
| M7 | Both axes from {1, R1} | Variant with sensitivity Level 1 + resistance R1 | Sensitive=Level1, Resistance=LevelR1, Combined=LevelR1 | Source 4 + SOP separate axes |
| M8 | No associations → NotActionable | Variant with empty association list | Combined/Sensitive/Resistance all None; IsActionable=false | Source 4 + ASSUMPTION A1 |
| M9 | Single Level 1 | Variant with one Level-1 association | Combined=Level1, Sensitive=Level1, Resistance=None, IsActionable=true | Source 2/4 |
| M10 | AssessActionability preserves order & count | List of 3 variants | 3 assessments, same order | INV-5 |
| M11 | ClassifyActionabilityLevel single highest | {3B, 2, 4} on one variant | Returns Level2 (highest combined) | Source 4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null variant list | `AssessActionability(null)` | `ArgumentNullException` | Library convention |
| S2 | Null associations on a variant | construct input with null association list | `ArgumentNullException` at construction/assess | Defensive validation |
| S3 | GetTherapyRecommendations ordering | variant with {4, 1, 3A} drug rows | rows ordered Level1, Level3A, Level4 (desc) | Source 4 order |
| S4 | GetTherapyRecommendations empty | variant with no associations | empty list (not null) | Annotator empty HIGHEST_LEVEL |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Standard-care grouping | IsStandardCare for levels 1,2,R1 true; 3A,3B,4,R2 false | per SOP grouping | SOP v3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing implementation of `AssessActionability` / OncoKB levels in `OncologyAnalyzer.cs` (only AMP/ASCO/CAP tiering for ONCO-ANNOT-001).
- No existing test file `OncologyAnalyzer_AssessActionability_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| M8 | ❌ Missing | new unit |
| M9 | ❌ Missing | new unit |
| M10 | ❌ Missing | new unit |
| M11 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AssessActionability_Tests.cs` — all cases above.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_AssessActionability_Tests.cs` | Canonical for ONCO-ACTION-001 | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | S1 | ❌ Missing | Implemented | ✅ Done |
| 13 | S2 | ❌ Missing | Implemented | ✅ Done |
| 14 | S3 | ❌ Missing | Implemented | ✅ Done |
| 15 | S4 | ❌ Missing | Implemented | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `ClassifyActionabilityLevel_CombinedOrder_*` |
| M2 | ✅ Covered | `ClassifyActionabilityLevel_R1AndLevel1_ReturnsR1` |
| M3 | ✅ Covered | `ClassifyActionabilityLevel_Level4AndR2_ReturnsLevel4` |
| M4 | ✅ Covered | `AssessActionability_SensitiveLevel2And3A_ReturnsLevel2` |
| M5 | ✅ Covered | `AssessActionability_Sensitive3A3B4_ReturnsLevel3A` |
| M6 | ✅ Covered | `AssessActionability_ResistanceR1R2_ReturnsR1` |
| M7 | ✅ Covered | `AssessActionability_BothAxes_ReportsEachAndCombined` |
| M8 | ✅ Covered | `AssessActionability_NoAssociations_NotActionable` |
| M9 | ✅ Covered | `AssessActionability_SingleLevel1_*` |
| M10 | ✅ Covered | `AssessActionability_PreservesOrderAndCount` |
| M11 | ✅ Covered | `ClassifyActionabilityLevel_MixedSet_ReturnsHighest` |
| S1 | ✅ Covered | `AssessActionability_NullVariants_Throws` |
| S2 | ✅ Covered | `VariantActionabilityInput_NullAssociations_Throws` |
| S3 | ✅ Covered | `GetTherapyRecommendations_OrdersByDescendingLevel` |
| S4 | ✅ Covered | `GetTherapyRecommendations_NoAssociations_ReturnsEmpty` |
| C1 | ✅ Covered | `IsStandardCare_LevelGrouping_MatchesSop` |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | No leveled association → `NotActionable` (no level); the level enum value name is ours, behavior matches annotator empty HIGHEST_LEVEL. | M8, S4, INV-4 |
| A2 | Drug–gene–level knowledgebase is caller-supplied; library ranks levels only (framework boundary). | all |

---

## 7. Open Questions / Decisions

1. Decision: scoped to OncoKB therapeutic levels (Chakravarty 2017) because ONCO-ANNOT-001 already covers AMP/ASCO/CAP tiering; no overlap.
2. Decision: implemented in `OncologyAnalyzer` (per session scope) rather than a new `ClinicalInterpreter` class named in the Registry; Registry method index updated to reflect actual class.
