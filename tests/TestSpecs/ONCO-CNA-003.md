# Test Specification: ONCO-CNA-003

**Test Unit ID:** ONCO-CNA-003
**Area:** Oncology
**Algorithm:** Homozygous (Deep) Deletion Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Cheng et al. (2017) Pan-cancer homozygous deletions, Nat Commun 8:1221 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5663922/ | 2026-06-14 |
| 2 | cBioPortal — Discrete Copy Number file format | 5 | https://docs.cbioportal.org/file-formats/ | 2026-06-14 |
| 3 | cBioPortal — FAQ (Deep/Shallow Deletion, −2..2) | 5 | https://docs.cbioportal.org/user-guide/faq/ | 2026-06-14 |
| 4 | CNVkit `cnvlib/call.py` `absolute_threshold` (integer CN; via ONCO-CNA-001) | 3 | https://cnvkit.readthedocs.io/ | 2026-06-14 |
| 5 | NCBI Gene — TP53/RB1/CDKN2A/PTEN/BRCA1/BRCA2 cytogenetic locations | 5 | https://www.ncbi.nlm.nih.gov/gene/7157 (+5) | 2026-06-14 |

### 1.2 Key Evidence Points

1. A homozygous deletion is a region with **zero copies of both alleles** → total/absolute copy number 0 — Cheng et al. (2017).
2. cBioPortal discrete scale: **−2 = Deep Deletion = (possible) homozygous deletion**; −1 = shallow / heterozygous loss (NOT homozygous); 0 diploid; 1 gain; 2 amplification — cBioPortal File-Formats / FAQ.
3. The repository already maps **integer copy number 0 ⇒ `CopyNumberState.DeepDeletion`** (log2 ≤ −1.1 by CNVkit `absolute_threshold` default) — ONCO-CNA-001 / CNVkit.
4. Tumour-suppressor arms (NCBI Gene): TP53 17p, RB1 13q, CDKN2A 9p, PTEN 10q, BRCA1 17q, BRCA2 13q.

### 1.3 Documented Corner Cases

- A single-copy loss (−1 / integer CN 1) is heterozygous, NOT a homozygous deletion (cBioPortal; Cheng et al. — one allele remains).
- Boundary: CNVkit assigns CN by "less than or equal to each threshold in sequence", so log2 exactly at the deletion cutoff (−1.1) is CN 0 (homozygous).

### 1.4 Known Failure Modes / Pitfalls

1. Reporting a shallow/heterozygous loss as homozygous — source-distinguished by total CN 0 vs ≥1 (Cheng et al.).
2. Purity/ploidy effects make discrete calls putative — interpretation caveat, does not change the CN-0 definition (cBioPortal).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectHomozygousDeletions(segments, thresholds?, ploidy?)` | OncologyAnalyzer | **Canonical** | Reports segments whose integer CN = 0 (DeepDeletion), order-preserving filter |
| `IsHomozygousDeletion(segment, thresholds?, ploidy?)` | OncologyAnalyzer | **Internal** | Predicate: integer CN == 0 |
| `IdentifyDeletedTumorSuppressors(deletions)` | OncologyAnalyzer | **Canonical** | Arm → gene panel (NCBI Gene) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A segment is reported iff its classified integer copy number is exactly 0 (DeepDeletion) | Yes | Cheng et al. (CN 0); cBioPortal (−2); CNVkit integer CN |
| INV-2 | A single-copy loss (integer CN 1) is never reported as homozygous | Yes | cBioPortal (−1 ≠ −2); Cheng et al. |
| INV-3 | Result is a subset of the input in input order (order/idempotence-preserving filter) | Yes | Filter semantics (mirror DetectFocalAmplifications) |
| INV-4 | Tumour-suppressor mapping is by chromosome arm, panel order, each gene once | Yes | NCBI Gene arms |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Deep deletion detected | Segment log2 = −2.0 (default thresholds → CN 0) | Reported as homozygous deletion | Cheng et al. (CN 0); cBioPortal −2; CNVkit |
| M2 | Single-copy loss excluded | Segment log2 = −0.5 (CN 1) | NOT reported | cBioPortal −1; Cheng et al. |
| M3 | Neutral excluded | Segment log2 = 0.0 (CN 2) | NOT reported | cBioPortal 0 |
| M4 | Gain/amplification excluded | log2 = 0.5 (CN 3) and log2 = 1.0 (CN ≥4) | NOT reported | cBioPortal 1 / 2 |
| M5 | Mixed set filtered in order | Segments [CN0, CN1, CN0, CN2] | Only the two CN-0 segments, in input order | INV-1, INV-3 |
| M6 | Boundary at deletion cutoff | log2 = −1.1 exactly (≤ −1.1 → CN 0) | Reported (homozygous) | CNVkit "≤ each threshold" |
| M7 | Just above cutoff not homozygous | log2 = −1.0999 (> −1.1 → CN 1) | NOT reported | CNVkit threshold |
| M8 | IsHomozygousDeletion predicate | CN-0 segment true, CN-1 segment false | true / false | INV-1, INV-2 |
| M9 | Map TP53 (17p) | Homozygous deletion on 17p | {"TP53"} | NCBI Gene TP53 17p13.1 |
| M10 | Map RB1 + BRCA2 (13q) | Homozygous deletion on 13q | {"RB1","BRCA2"} (panel order) | NCBI Gene RB1 13q14.2, BRCA2 13q13.1 |
| M11 | Map CDKN2A (9p) | Homozygous deletion on 9p | {"CDKN2A"} | NCBI Gene CDKN2A 9p21.3 |
| M12 | Map PTEN (10q) | Homozygous deletion on 10q | {"PTEN"} | NCBI Gene PTEN 10q23.31 |
| M13 | Map BRCA1 (17q) | Homozygous deletion on 17q | {"BRCA1"} | NCBI Gene BRCA1 17q21.31 |
| M14 | Custom thresholds shift CN-0 | Raise deletion cutoff so a log2=−0.5 segment becomes CN 0 | Reported under custom thresholds | CNVkit thresholds are parameters |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Non-panel arm yields no gene | Homozygous deletion on an arm with no panel gene (e.g. 1p) | empty | gene panel is closed |
| S2 | Distinct genes once | Two 13q deletions | {"RB1","BRCA2"} not duplicated | INV-4 |
| S3 | Ploidy parameter | Triploid reference shifts CN-0 boundary | classification respects ploidy | CNVkit n = ploidy·2^log2 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty input | empty segments | empty result | mirror sibling |
| C2 | Null segments | null | ArgumentNullException | mirror sibling |
| C3 | Invalid segment | ArmLength ≤ 0 or End ≤ Start | ArgumentException | mirror ValidateArmSegment |
| C4 | Null deletions to mapper | null | ArgumentNullException | mirror IdentifyAmplifiedOncogenes |
| C5 | NaN log2 no-call | log2 = NaN | NOT reported (neutral no-call, CN = ploidy) | CNVkit NaN → neutral |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `DetectHomozygousDeletions` / `IdentifyDeletedTumorSuppressors` (grep of OncologyAnalyzer + tests directory). Brand-new unit.
- Sibling pattern: `OncologyAnalyzer_DetectFocalAmplifications_Tests.cs` (ONCO-CNA-002).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14, S1–S3, C1–C5 | ❌ Missing | Brand-new unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectHomozygousDeletions_Tests.cs` — all cases for both canonical methods + the internal predicate.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_DetectHomozygousDeletions_Tests.cs | Canonical (this unit) | 22 |

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
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | S1 | ❌ Missing | Implemented | ✅ Done |
| 16 | S2 | ❌ Missing | Implemented | ✅ Done |
| 17 | S3 | ❌ Missing | Implemented | ✅ Done |
| 18 | C1 | ❌ Missing | Implemented | ✅ Done |
| 19 | C2 | ❌ Missing | Implemented | ✅ Done |
| 20 | C3 | ❌ Missing | Implemented | ✅ Done |
| 21 | C4 | ❌ Missing | Implemented | ✅ Done |
| 22 | C5 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 22
**✅ Done:** 22 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | DetectHomozygousDeletions_DeepDeletionSegment_IsReported |
| M2 | ✅ Covered | DetectHomozygousDeletions_SingleCopyLoss_NotReported |
| M3 | ✅ Covered | DetectHomozygousDeletions_NeutralSegment_NotReported |
| M4 | ✅ Covered | DetectHomozygousDeletions_GainAndAmplification_NotReported |
| M5 | ✅ Covered | DetectHomozygousDeletions_MixedSet_ReturnsCn0InOrder |
| M6 | ✅ Covered | DetectHomozygousDeletions_Log2AtDeletionCutoff_IsReported |
| M7 | ✅ Covered | DetectHomozygousDeletions_Log2JustAboveCutoff_NotReported |
| M8 | ✅ Covered | IsHomozygousDeletion_Cn0AndCn1_TrueThenFalse |
| M9 | ✅ Covered | IdentifyDeletedTumorSuppressors_Arm17p_MapsTp53 |
| M10 | ✅ Covered | IdentifyDeletedTumorSuppressors_Arm13q_MapsRb1AndBrca2 |
| M11 | ✅ Covered | IdentifyDeletedTumorSuppressors_Arm9p_MapsCdkn2a |
| M12 | ✅ Covered | IdentifyDeletedTumorSuppressors_Arm10q_MapsPten |
| M13 | ✅ Covered | IdentifyDeletedTumorSuppressors_Arm17q_MapsBrca1 |
| M14 | ✅ Covered | DetectHomozygousDeletions_CustomThresholds_ShiftsCn0Boundary |
| S1 | ✅ Covered | IdentifyDeletedTumorSuppressors_NonPanelArm_ReturnsEmpty |
| S2 | ✅ Covered | IdentifyDeletedTumorSuppressors_DuplicateArm_GenesOnce |
| S3 | ✅ Covered | DetectHomozygousDeletions_TriploidPloidy_RespectsBoundary |
| C1 | ✅ Covered | DetectHomozygousDeletions_EmptyInput_ReturnsEmpty |
| C2 | ✅ Covered | DetectHomozygousDeletions_NullInput_Throws |
| C3 | ✅ Covered | DetectHomozygousDeletions_InvalidSegment_Throws |
| C4 | ✅ Covered | IdentifyDeletedTumorSuppressors_NullInput_Throws |
| C5 | ✅ Covered | DetectHomozygousDeletions_NaNLog2_NotReported |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Tumour-suppressor panel (TP53/RB1/CDKN2A/PTEN/BRCA1/BRCA2) is a registry-supplied curated list; arms are source-backed (NCBI Gene) | IdentifyDeletedTumorSuppressors |
| 2 | Homozygous deletion = integer CN 0 via existing DeepDeletion classification (no new threshold invented) | DetectHomozygousDeletions, IsHomozygousDeletion |

---

## 7. Open Questions / Decisions

1. **Class placement / naming conflict (carried from ONCO-CNA-001/002):** registry names `CopyNumberAnalyzer`; the implemented analyzer is `OncologyAnalyzer` (the `CopyNumberAnalyzer` names were superseded in ONCO-CNA-001). Methods are added to `OncologyAnalyzer` for consistency. Recorded here per the checklist-conflict rule.
2. **Segment type reuse:** reuses `CopyNumberArmSegment` (introduced in ONCO-CNA-002, carries Arm + ArmLength + Log2Ratio) so the same input drives both focal-amplification and homozygous-deletion detection. No new segment type needed.
