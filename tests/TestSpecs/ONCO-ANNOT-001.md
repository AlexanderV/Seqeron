# Test Specification: ONCO-ANNOT-001

**Test Unit ID:** ONCO-ANNOT-001
**Area:** Oncology
**Algorithm:** Cancer-Specific Variant Annotation (AMP/ASCO/CAP 2017 four-tier clinical-significance classification)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Li MM et al. (2017). Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer (AMP/ASCO/CAP). J Mol Diagn 19(1):4–23. | 1 | https://doi.org/10.1016/j.jmoldx.2016.10.002 (full text read from https://ocpe.mcw.edu/sites/default/files/course/2024-03/AMP-ASCO-CAP%20guidelines%20-%20somatic%20variants.pdf) | 2026-06-14 |
| 2 | Tate JG et al. (2019). COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res 47(D1):D941–D947. | 1/5 | https://doi.org/10.1093/nar/gky1015 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Four tiers: I (strong clinical significance, Level A/B), II (potential, Level C/D), III (unknown), IV (benign / likely benign) — Li et al. (2017) Figure 2.
2. Tier I = Level A and B evidence; Tier II = Level C or D evidence (Figure 2, verbatim).
3. Benign population cutoff: "the work group recommends using 1% (0.01) as a primary cutoff" for eliminating polymorphic/benign variants — Li et al. (2017), Population Databases.
4. Tier IV population criterion: "MAF ≥ 1% in the general population; or high MAF in some ethnic populations" — Table 7. Tier IV also: "No existing published evidence of cancer association" (Figure 2).
5. Tier III: rare (absent/extremely low MAF), no clinical evidence, but no convincing cancer association yet not common — Table 6 / Figure 2.
6. COSMIC is an external curated somatic-mutation database; lookups must be against caller-supplied records — Tate et al. (2019).

### 1.3 Documented Corner Cases

- MAF ≥ 1% with no clinical evidence ⇒ Tier IV (Table 7).
- Significant biomarker (Level A/B) that also appears in population databases ⇒ still Tier I (assigned by evidence level, Figure 2).
- Rare variant with cancer association but no clinical evidence ⇒ Tier III; if no cancer association ⇒ Tier IV (Figure 2 boxes).
- COSMIC catalog miss ⇒ return null, do not fabricate.

### 1.4 Known Failure Modes / Pitfalls

1. Downgrading a Level A/B biomarker to benign because it is common — wrong; the guideline categorizes by evidence level — Li et al. (2017) Figure 2.
2. Treating the 1% cutoff as strict `>` rather than `≥` — Table 7 reads "MAF ≥ 1%" / the 1% primary cutoff is inclusive for benign.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifyVariantTier(CancerVariantAnnotationInput)` | OncologyAnalyzer | Canonical | Core Figure 2 decision rule (per-variant). |
| `AnnotateCancerVariants(IEnumerable<CancerVariantAnnotationInput>)` | OncologyAnalyzer | Canonical | Batch wrapper; one annotation per variant, input order. |
| `GetCOSMICAnnotation(CancerVariantAnnotationInput, IReadOnlyDictionary<(string,string),string>)` | OncologyAnalyzer | Canonical | Exact-match lookup against caller-supplied COSMIC catalog. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every variant maps to exactly one of the four tiers (total, deterministic). | Yes | Li et al. (2017) Figure 2 (four exhaustive categories) |
| INV-2 | Level A/B ⇒ Tier I; Level C/D ⇒ Tier II, regardless of MAF / cancer association. | Yes | Li et al. (2017) Figure 2 |
| INV-3 | With no evidence level: MAF ≥ 0.01 OR no cancer association ⇒ Tier IV; else Tier III. | Yes | Li et al. (2017) Tables 6/7, Figure 2 |
| INV-4 | `AnnotateCancerVariants` output count = input count, in input order. | Yes | composition of per-variant rule |
| INV-5 | `GetCOSMICAnnotation` returns a value iff (gene, protein change) is in the supplied catalog. | Yes | Tate et al. (2019); external lookup |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Level A ⇒ Tier I | Variant with Level A evidence | `TierI_StrongClinicalSignificance` | Li (2017) Fig 2 |
| M2 | Level B ⇒ Tier I | Variant with Level B evidence | `TierI_StrongClinicalSignificance` | Li (2017) Fig 2 |
| M3 | Level C ⇒ Tier II | Variant with Level C evidence | `TierII_PotentialClinicalSignificance` | Li (2017) Fig 2 |
| M4 | Level D ⇒ Tier II | Variant with Level D evidence | `TierII_PotentialClinicalSignificance` | Li (2017) Fig 2 |
| M5 | Common SNP ⇒ Tier IV | No evidence, MAF 0.25 | `TierIV_BenignOrLikelyBenign` | Li (2017) Table 7 |
| M6 | No evidence, low MAF, no cancer assoc. ⇒ Tier IV | MAF 0.0001, association false | `TierIV_BenignOrLikelyBenign` | Li (2017) Fig 2 Tier IV box |
| M7 | Rare VUS ⇒ Tier III | No evidence, MAF 0.0001, association true | `TierIII_UnknownClinicalSignificance` | Li (2017) Table 6 / Fig 2 |
| M8 | Evidence priority over MAF | Level A, MAF 0.30, association false | `TierI_StrongClinicalSignificance` | Li (2017) Fig 2 (by evidence level) |
| M9 | MAF boundary = 0.01 ⇒ Tier IV | No evidence, MAF exactly 0.01, association true | `TierIV_BenignOrLikelyBenign` | Li (2017) "≥ 1%" cutoff |
| M10 | MAF just below 0.01 ⇒ Tier III | No evidence, MAF 0.0099, association true | `TierIII_UnknownClinicalSignificance` | Li (2017) "≥ 1%" cutoff (inclusive) |
| M11 | `AnnotateCancerVariants` batch | Mixed batch (A, C, common SNP, rare VUS) | Tiers [I, II, IV, III], same order, count = 4 | Li (2017); INV-4 |
| M12 | COSMIC hit | Lookup (BRAF, p.V600E) in catalog | returns catalog value "COSV56056643" | Tate (2019) |
| M13 | COSMIC miss | Lookup (TP53, p.R175H) not in catalog | returns null | Tate (2019) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Negative MAF | MAF = -0.1 | `ArgumentOutOfRangeException` | input validation |
| S2 | MAF > 1 | MAF = 1.5 | `ArgumentOutOfRangeException` | input validation |
| S3 | NaN MAF | MAF = NaN | `ArgumentOutOfRangeException` | input validation |
| S4 | Null variants | `AnnotateCancerVariants(null)` | `ArgumentNullException` | API contract |
| S5 | Null catalog | `GetCOSMICAnnotation(v, null)` | `ArgumentNullException` | API contract |
| S6 | Empty batch | `AnnotateCancerVariants([])` | empty list | boundary |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Totality | Every (level × maf-band × assoc) combination yields a defined tier | no exceptions; INV-1 holds | property check |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests exist for `AnnotateCancerVariants`, `ClassifyVariantTier`, or `GetCOSMICAnnotation` (these methods are newly added in this unit). Existing `OncologyAnalyzer_*` test files cover unrelated methods (somatic calling, VAF, drivers, artifacts).
- Canonical file created: `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnnotateCancerVariants_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M13 | ❌ Missing | New methods; no prior tests. |
| S1–S6 | ❌ Missing | New methods; no prior tests. |
| C1 | ❌ Missing | New property test. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnnotateCancerVariants_Tests.cs` — all ONCO-ANNOT-001 tests.
- **Remove:** none (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_AnnotateCancerVariants_Tests.cs` | Canonical (ONCO-ANNOT-001) | 20 |

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
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | S1 | ❌ Missing | Implemented | ✅ Done |
| 15 | S2 | ❌ Missing | Implemented | ✅ Done |
| 16 | S3 | ❌ Missing | Implemented | ✅ Done |
| 17 | S4 | ❌ Missing | Implemented | ✅ Done |
| 18 | S5 | ❌ Missing | Implemented | ✅ Done |
| 19 | S6 | ❌ Missing | Implemented | ✅ Done |
| 20 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 20
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `ClassifyVariantTier_LevelA_ReturnsTierI` |
| M2 | ✅ Covered | `ClassifyVariantTier_LevelB_ReturnsTierI` |
| M3 | ✅ Covered | `ClassifyVariantTier_LevelC_ReturnsTierII` |
| M4 | ✅ Covered | `ClassifyVariantTier_LevelD_ReturnsTierII` |
| M5 | ✅ Covered | `ClassifyVariantTier_CommonPolymorphism_ReturnsTierIV` |
| M6 | ✅ Covered | `ClassifyVariantTier_RareNoCancerAssociation_ReturnsTierIV` |
| M7 | ✅ Covered | `ClassifyVariantTier_RareWithCancerAssociation_ReturnsTierIII` |
| M8 | ✅ Covered | `ClassifyVariantTier_LevelAWithHighMaf_StillReturnsTierI` |
| M9 | ✅ Covered | `ClassifyVariantTier_MafExactlyAtOnePercent_ReturnsTierIV` |
| M10 | ✅ Covered | `ClassifyVariantTier_MafJustBelowOnePercent_ReturnsTierIII` |
| M11 | ✅ Covered | `AnnotateCancerVariants_MixedBatch_PreservesOrderAndTiers` |
| M12 | ✅ Covered | `GetCOSMICAnnotation_VariantInCatalog_ReturnsCosmicId` |
| M13 | ✅ Covered | `GetCOSMICAnnotation_VariantNotInCatalog_ReturnsNull` |
| S1 | ✅ Covered | `ClassifyVariantTier_NegativeMaf_Throws` |
| S2 | ✅ Covered | `ClassifyVariantTier_MafAboveOne_Throws` |
| S3 | ✅ Covered | `ClassifyVariantTier_NaNMaf_Throws` |
| S4 | ✅ Covered | `AnnotateCancerVariants_NullVariants_Throws` |
| S5 | ✅ Covered | `GetCOSMICAnnotation_NullCatalog_Throws` |
| S6 | ✅ Covered | `AnnotateCancerVariants_EmptyBatch_ReturnsEmpty` |
| C1 | ✅ Covered | `ClassifyVariantTier_AllEvidenceAndMafCombinations_ProduceDefinedTier` |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Evidence inputs (level, MAF, cancer association) are caller-supplied; the library does not reproduce curated databases. | Method contract; all tests |
| 2 | Tier III vs IV discriminator: no evidence level + (MAF ≥ 0.01 OR no cancer association) ⇒ Tier IV, else Tier III (direct reading of Figure 2 boxes / Tables 6–7). | M6, M7, M9, M10, INV-3 |

---

## 7. Open Questions / Decisions

1. **COSMIC lookup shape:** COSMIC content is external and cannot be hardcoded (Tate 2019); the lookup is against a caller-supplied dictionary keyed by (gene, protein change). Decided — consistent with the existing `MatchCancerHotspots` caller-supplied-set pattern in this class.
2. **No other open questions.** Both assumptions are direct readings of the cited guideline, not unresolved correctness gaps.
