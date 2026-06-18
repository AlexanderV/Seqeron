# Test Specification: ONCO-CNA-002

**Test Unit ID:** ONCO-CNA-002
**Area:** Oncology
**Algorithm:** Focal Amplification Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mermel et al. (2011) GISTIC2.0, Genome Biology 12:R41 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/ | 2026-06-14 |
| 2 | Broad Institute GISTIC2 docs (`broad_len_cutoff`, `t_amp`) | 3 | https://broadinstitute.github.io/gistic2/ | 2026-06-14 |
| 3 | CNVkit — Calling copy number gains and losses | 3 | https://cnvkit.readthedocs.io/en/stable/calling.html | 2026-06-14 |
| 4 | NCBI Gene (ERBB2 2064, MYC 4609, EGFR 1956, CCND1 595, MDM2 4193, CDK4 1019) | 5 | https://www.ncbi.nlm.nih.gov/gene/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. A focal SCNA has length < 98% of a chromosome arm; an event occupying > 98% of an arm is arm-level/broad. — Mermel et al. (2011).
2. The cutoff is the configurable `broad_len_cutoff`, default 0.98 (fraction of chromosome arm). — GISTIC2 docs.
3. A segment is "amplified" when its copy-number gain exceeds `t_amp`, default 0.1 (log2). — GISTIC2 docs.
4. A single-copy gain is log2(3/2) = 0.585, well above the 0.1 amplitude threshold. — CNVkit docs.
5. Oncogene arms: ERBB2 17q, MYC 8q, EGFR 7p, CCND1 11q, MDM2 12q, CDK4 12q. — NCBI Gene.

### 1.3 Documented Corner Cases

- A segment ≥ 98% of its arm is arm-level and must NOT be reported as focal even when highly amplified (Mermel et al. 2011).
- A gain not exceeding `t_amp` (0.1) is not amplified and is excluded regardless of length (GISTIC2 docs).

### 1.4 Known Failure Modes / Pitfalls

1. Treating a whole-arm amplification as focal (off-by-one on the 98% rule). — Mermel et al. (2011).
2. Reporting low-level (artifactual) amplitude segments as amplifications. — GISTIC2 `t_amp`.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectFocalAmplifications(segments, thresholds?)` | OncologyAnalyzer | Canonical | GISTIC2 length(<0.98·arm) + amplitude(>t_amp) predicate |
| `IdentifyAmplifiedOncogenes(amplifications)` | OncologyAnalyzer | Canonical | Maps focal amplifications to ERBB2/MYC/EGFR/CCND1/MDM2/CDK4 by chromosome arm |
| `IsFocalAmplification(segment, thresholds)` | OncologyAnalyzer | Internal | Single-segment predicate; tested via canonicals |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported focal amplification has length/armLength < 0.98 (strictly below cutoff) | Yes | Mermel et al. (2011); GISTIC2 `broad_len_cutoff` 0.98 |
| INV-2 | Every reported focal amplification has log2 > t_amp (default 0.1) | Yes | GISTIC2 `t_amp` 0.1 |
| INV-3 | DetectFocalAmplifications preserves a subset of the input segments (no fabricated segments); order = input order | Yes | Filtering semantics |
| INV-4 | An oncogene is reported only for an arm that carries a focal amplification | Yes | Mapping operates on focal amplifications only |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Focal high-amp segment | 17q, arm 1e6, seg 5e5 (0.50), log2 1.0 | Reported as focal amplification | Mermel et al. (2011) <0.98; `t_amp` 0.1 |
| M2 | Whole-arm amplification | 8q, arm 1e6, seg 9.9e5 (0.99), log2 1.5 | NOT focal (arm-level) | Mermel et al. (2011) >98% ⇒ arm-level |
| M3 | Boundary exactly 0.98 | 11q, arm 1e6, seg 9.8e5 (0.98), log2 1.0 | NOT focal (not strictly < 0.98) | GISTIC2 `broad_len_cutoff` 0.98 |
| M4 | Low-amplitude focal segment | 7p, arm 1e6, seg 3e5 (0.30), log2 0.05 | NOT reported (≤ t_amp) | GISTIC2 `t_amp` 0.1 |
| M5 | Just above amplitude cutoff | arm 1e6, seg 1e5 (0.10), log2 0.585 (single-copy gain) | Reported as focal amplification | CNVkit log2(3/2)=0.585 > 0.1 |
| M6 | Map 17q → ERBB2 | focal amp on 17q | Oncogene set contains ERBB2 | NCBI Gene 2064 (17q12) |
| M7 | Map 8q → MYC | focal amp on 8q | Oncogene set contains MYC | NCBI Gene 4609 (8q24.21) |
| M8 | Map 7p → EGFR | focal amp on 7p | Oncogene set contains EGFR | NCBI Gene 1956 (7p11.2) |
| M9 | Map 11q → CCND1 | focal amp on 11q | Oncogene set contains CCND1 | NCBI Gene 595 (11q13.3) |
| M10 | Map 12q → MDM2 and CDK4 | focal amp on 12q | Oncogene set contains both MDM2 and CDK4 | NCBI Gene 4193 (12q15), 1019 (12q14.1) |
| M11 | Order/subset preservation | mixed list (focal, arm-level, low-amp, focal) | Output = the two focal amps in input order | INV-3 |
| M12 | Mapping ignores non-amplified arms | only low-amp segment on 17q | ERBB2 NOT reported | INV-4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Custom thresholds | t_amp 0.3, log2 0.2 segment | NOT reported (below custom t_amp) | Parameter override |
| S2 | Arm with no panel oncogene | focal amp on 5q | Oncogene set empty | Only panel arms map |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null segments | DetectFocalAmplifications(null) | ArgumentNullException | Guard |
| C2 | Empty segments | DetectFocalAmplifications([]) | empty list | Guard |
| C3 | Null amplifications | IdentifyAmplifiedOncogenes(null) | ArgumentNullException | Guard |
| C4 | Invalid arm length | armLength ≤ 0 | ArgumentException | Validation |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `DetectFocalAmplifications` / `IdentifyAmplifiedOncogenes`; grep over `src/` and `tests/` found no such methods. New unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M12, S1–S2, C1–C4 | ❌ Missing | New unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectFocalAmplifications_Tests.cs` — all cases for both canonical methods.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_DetectFocalAmplifications_Tests.cs` | Canonical | 18 |

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
| 13 | S1 | ❌ Missing | Implemented | ✅ Done |
| 14 | S2 | ❌ Missing | Implemented | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |
| 17 | C3 | ❌ Missing | Implemented | ✅ Done |
| 18 | C4 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 18
**✅ Done:** 18 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | DetectFocalAmplifications_FocalHighAmp_Reported |
| M2 | ✅ Covered | DetectFocalAmplifications_WholeArm_NotReported |
| M3 | ✅ Covered | DetectFocalAmplifications_ExactlyCutoff_NotReported |
| M4 | ✅ Covered | DetectFocalAmplifications_LowAmplitude_NotReported |
| M5 | ✅ Covered | DetectFocalAmplifications_JustAboveAmpCutoff_Reported |
| M6 | ✅ Covered | IdentifyAmplifiedOncogenes_Arm17q_ReturnsErbb2 |
| M7 | ✅ Covered | IdentifyAmplifiedOncogenes_Arm8q_ReturnsMyc |
| M8 | ✅ Covered | IdentifyAmplifiedOncogenes_Arm7p_ReturnsEgfr |
| M9 | ✅ Covered | IdentifyAmplifiedOncogenes_Arm11q_ReturnsCcnd1 |
| M10 | ✅ Covered | IdentifyAmplifiedOncogenes_Arm12q_ReturnsMdm2AndCdk4 |
| M11 | ✅ Covered | DetectFocalAmplifications_MixedList_PreservesFocalSubsetInOrder |
| M12 | ✅ Covered | IdentifyAmplifiedOncogenes_NonAmplifiedArm_NotMapped (via DetectFocal empty) |
| S1 | ✅ Covered | DetectFocalAmplifications_CustomTamp_BelowThreshold_NotReported |
| S2 | ✅ Covered | IdentifyAmplifiedOncogenes_ArmWithoutPanelGene_Empty |
| C1 | ✅ Covered | DetectFocalAmplifications_Null_Throws |
| C2 | ✅ Covered | DetectFocalAmplifications_Empty_ReturnsEmpty |
| C3 | ✅ Covered | IdentifyAmplifiedOncogenes_Null_Throws |
| C4 | ✅ Covered | DetectFocalAmplifications_NonPositiveArmLength_Throws |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Amplitude "amplified" test uses GISTIC2 `t_amp` = 0.1 combined with the paper's length rule | DetectFocalAmplifications predicate |
| 2 | Arm label + arm length supplied by caller (no bundled cytoband table) | Segment input contract |

---

## 7. Open Questions / Decisions

1. Registry names the class `CopyNumberAnalyzer`; per ONCO-CNA-001 the `CopyNumberAnalyzer`/`CallCopyNumberStates` names were superseded by `OncologyAnalyzer`. This unit follows the established `OncologyAnalyzer` placement (consistent with ONCO-CNA-001). Conflict noted; registry method-index rows point to ONCO-CNA-002 already.
