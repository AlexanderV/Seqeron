# TestSpec: CRISPR-GUIDE-001 - Guide RNA Design

## Test Unit Identification
- **Test Unit ID**: CRISPR-GUIDE-001
- **Algorithm Group**: MolTools
- **Algorithm Name**: Guide_RNA_Design
- **Canonical Methods**:
  - `CrisprDesigner.DesignGuideRnas(DnaSequence, int, int, CrisprSystemType, GuideRnaParameters?)`
  - `CrisprDesigner.EvaluateGuideRna(string, CrisprSystemType, GuideRnaParameters?)`
- **Related Documentation**: [Guide_RNA_Design.md](../../docs/algorithms/MolTools/Guide_RNA_Design.md)

## Evidence Summary

### Primary Sources
| Source | Key Information | URL |
|--------|-----------------|-----|
| Addgene CRISPR Guide | gRNA structure, 20bp spacer, seed sequence (8-10bp at 3'), PAM requirement, Pol III TTTT terminator | https://www.addgene.org/guides/crispr/ |
| Wikipedia: Guide RNA | GC >50% optimal, length 17-24bp standard 20bp | https://en.wikipedia.org/wiki/Guide_RNA |
| Wikipedia: PAM | SpCas9 PAM is 5'-NGG-3', PAM required for cleavage | https://en.wikipedia.org/wiki/Protospacer_adjacent_motif |

### Evidence-Backed Parameters
| Parameter | Evidence Value | Source | Implementation Value | Status |
|-----------|----------------|--------|---------------------|--------|
| Guide length | 20bp standard (17-24bp range) | Wikipedia Guide RNA | 20bp | ✅ Aligned |
| Optimal GC | >50% | Wikipedia Guide RNA | 40-70% range | ✅ Broader range acceptable per common tools |
| Seed region | 8-10bp at 3' end | Addgene CRISPR Guide | 10bp (upper bound) | ✅ Aligned |
| Poly-T termination | TTTT (Pol III terminator) | Addgene CRISPR Guide | TTTT detection | ✅ Aligned |
| SpCas9 PAM | 5'-NGG-3' | Addgene, Wikipedia PAM | NGG | ✅ Aligned |
| SaCas9 PAM | NNGRRT | Addgene PAM table | NNGRRT | ✅ Aligned |
| Cas12a PAM | 5'-TTTV-3' | Addgene PAM table | TTTV | ✅ Aligned |
| Scaffold | tracrRNA-derived scaffold | Addgene: sgRNA = spacer + scaffold | SpCas9 76nt scaffold | ✅ Aligned |

### Scoring Model (verified via exact-value tests)
| Penalty | Formula | Verified By |
|---------|---------|-------------|
| GC outside [MinGc, MaxGc] | −(delta × 2) per % | M-001 (0 penalty), M-002 (−80), M-003 (−60) |
| Poly-T (TTTT+) | −20 flat | M-004 (score 80), S-008 (score 80) |
| Self-complementarity > 0.3 | −(selfComp × 30) | C-001 (−9.375 on 8bp palindrome) |
| Seed GC outside [30%, 80%] | −5 flat | S-010 (−5), S-011 (−5), M-002 (combined), M-003 (combined) |
| Restriction site (common) | −5 flat | S-001 (score 95) |
| Score floor | max(0, score) | C-002 (clamped to 0) |

---

## Test Categories

### MUST Tests (Critical Functionality) — 9 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| M-001 | `EvaluateGuideRna_OptimalGuide_HighScore` | "ACGTACGTACGTACGTACGT" (50% GC) | Score=100, GC=50, SeedGC=50, Issues=[], SelfComp=0.15 | ✅ |
| M-002 | `EvaluateGuideRna_LowGcContent_LowerScore` | "AAAAAAAAAAAAAAAAAAAA" (0% GC) | Score=15, GC=0, SeedGC=0, Issues=2 [Low GC, Suboptimal seed GC] | ✅ |
| M-003 | `EvaluateGuideRna_HighGcContent_LowerScore` | "GCGCGCGCGCGCGCGCGCGC" (100% GC) | Score=35, GC=100, SeedGC=100, Issues=2 [High GC, Suboptimal seed GC] | ✅ |
| M-004 | `EvaluateGuideRna_HasPolyT_Penalized` | "ACGTACGTTTTTACGTACGT" (TTTT) | Score=80, GC=40, HasPolyT=true, Issues=1 [TTTT] | ✅ |
| M-005 | `EvaluateGuideRna_EmptyGuide_ThrowsException` | "" | Throws ArgumentNullException | ✅ |
| M-006 | `EvaluateGuideRna_FullGuideRna_IncludesScaffold` | any 20bp guide | FullGuideRna = guide + scaffold (96 chars total) | ✅ |
| M-007 | `DesignGuideRnas_NullSequence_ThrowsException` | null | Throws ArgumentNullException | ✅ |
| M-008 | `DesignGuideRnas_InvalidRegionStart_ThrowsException` | regionStart = −1 | Throws ArgumentOutOfRangeException | ✅ |
| M-009 | `DesignGuideRnas_InvalidRegionEnd_ThrowsException` | regionEnd > len | Throws ArgumentOutOfRangeException | ✅ |

---

### SHOULD Tests (Important Quality) — 11 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| S-001 | `EvaluateGuideRna_RestrictionSite_Penalized` | "ACGTGAATTCACGTACGTAC" (EcoRI) | Score=95, GC=45, Issues=1 [restriction site] | ✅ |
| S-002 | `EvaluateGuideRna_CalculatesSeedGc` | "AAAAAAAAAAAAACGTACGT" | GC=20, SeedGC=40, Score=60, Issues=1 [Low GC] | ✅ |
| S-003 | `DesignGuideRnas_WithPamInRegion_ReturnsGuides` | seq with AGG PAM | Count=1, Pos=24, Score=100, Forward=true | ✅ |
| S-004 | `GuideRnaParameters_Default_HasValidValues` | Default | MinGC=40, MaxGC=70, MinScore=50, AvoidPolyT=true, CheckSelfComp=true | ✅ |
| S-005 | `GuideRnaParameters_CustomValues_Respected` | Custom(30,80,40,false,false) | All 5 fields preserved | ✅ |
| S-006 | `EvaluateGuideRna_BoundaryGc40Percent_NotPenalized` | "AAAAAAAAAAAAGCGCGCGC" (40%) | Score=100, GC=40, SeedGC=80, Issues=[] | ✅ |
| S-007 | `EvaluateGuideRna_BoundaryGc70Percent_NotPenalized` | "GCGCGCGCGCGCGCAAAAAA" (70%) | Score=100, GC=70, SeedGC=40, Issues=[] | ✅ |
| S-008 | `EvaluateGuideRna_ExactlyFourTs_TriggersPolyT` | "ACGTACGTACGATTTTACGT" (4 T's) | Score=80, HasPolyT=true, Issues=1 [TTTT] | ✅ |
| S-009 | `EvaluateGuideRna_ThreeConsecutiveTs_NoPolyT` | "ACGTACGTACGTACGTTTAC" (3 T's) | Score=100, HasPolyT=false, Issues=[] | ✅ |
| S-010 | `EvaluateGuideRna_SeedGcLow_Penalized` | "GCGCGCGCGCAAAAAAAAAA" | GC=50, SeedGC=0, Score=95, Issues=1 [Suboptimal seed GC] | ✅ |
| S-011 | `EvaluateGuideRna_SeedGcHigh_Penalized` | "AAAAAAAAAAGGGGGGGGGG" | GC=50, SeedGC=100, Score=95, Issues=1 [Suboptimal seed GC] | ✅ |

---

### COULD Tests (Edge Cases) — 10 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| C-001 | `EvaluateGuideRna_SelfComplementary_PenaltyTriggered` | "GCGCGCGC" (8bp, selfComp=0.3125) | Score=25.625, Issues ∋ "self-complementarity"; control "ACGTACGT" Score=100 | ✅ |
| C-002 | `EvaluateGuideRna_AllT_VeryLowScoreWithMultipleIssues` | "TTTTTTTTTTTTTTTTTTTT" | Score=0 (clamped), GC=0, SeedGC=0, HasPolyT=true, Issues=3 | ✅ |
| C-003 | `EvaluateGuideRna_NullGuide_ThrowsException` | null | Throws ArgumentNullException | ✅ |
| C-004 | `DesignGuideRnas_NoPamInRegion_ReturnsEmpty` | all-A sequence | Returns empty | ✅ |
| C-005 | `DesignGuideRnas_MultiplePams_ReturnsMultipleGuides` | 3 × (20bp+NGG) | Count=3, all Score=100 | ✅ |
| C-006 | `EvaluateGuideRna_SaCas9SystemType_ValidEvaluation` | "ACGTACGTACGTACGTACGT" SaCas9 | Score=100, System.Name="SaCas9", GuideLength=21 | ✅ |
| C-007 | `EvaluateGuideRna_BelowBoundaryGc_HasLowGcIssue` | "AAAAAAAAAAAAGCGCGCAT" (30%) | GC=30, Score=80, Issues=1 [Low GC] | ✅ |
| C-008 | `DesignGuideRnas_EntireSequenceAsRegion_Works` | 27bp seq with PAM | Count=1, Pos=4, Score=100 | ✅ |
| C-009 | `EvaluateGuideRna_AboveBoundaryGc_HasHighGcIssue` | "GCGCGCGCGCGCGCGCAAAA" (80%) | GC=80, Score=80, Issues=1 [High GC] | ✅ |
| C-010 | `DesignGuideRnas_MinScoreFiltering_ExcludesLowScoreGuides` | seq with PAM | MinScore=50 → Count=1; MinScore=101 → empty | ✅ |

---

## Property Tests (CrisprProperties.cs)

| Test | Assertion | Status |
|------|-----------|--------|
| `GuideRna_GcContent_InRange` | GC ∈ [0, 100] for all designed guides; test sequence produces guides | ✅ |
| `GuideRna_Score_InRange` | Score ∈ [0, 100] for all designed guides; test sequence produces guides | ✅ |

---

## Test File Mapping

| Test Category | Source File |
|--------------|-------------|
| Guide RNA Evaluation (30 tests) | CrisprDesigner_GuideRNA_Tests.cs |
| Guide RNA Properties (2 tests) | Properties/CrisprProperties.cs |
| Off-Target Analysis | CrisprDesigner_OffTarget_Tests.cs (CRISPR-OFF-001) |
| PAM Detection | CrisprDesigner_PAM_Tests.cs (CRISPR-PAM-001) |

---

## Coverage Summary

| Category | MUST | SHOULD | COULD | Total |
|----------|------|--------|-------|-------|
| Required | 9 | 11 | 10 | 30 |
| Implemented | 9 | 11 | 10 | 30 |
| Gap | 0 | 0 | 0 | **0** |

---

## Coverage Classification Audit (v3.0)

All code paths in `EvaluateGuideRna` and `DesignGuideRnas` are covered with exact-value assertions:

| Code Path | Covered By | Classification |
|-----------|-----------|---------------|
| GC penalty (low) | M-002 (Score=15) | ✅ Covered |
| GC penalty (high) | M-003 (Score=35) | ✅ Covered |
| GC boundary 40% (inclusive, no penalty) | S-006 (Score=100) | ✅ Covered |
| GC boundary 70% (inclusive, no penalty) | S-007 (Score=100) | ✅ Covered |
| GC below boundary 30% | C-007 (Score=80) | ✅ Covered |
| GC above boundary 80% | C-009 (Score=80) | ✅ Covered |
| Poly-T detection (4+ T's) | M-004, S-008 | ✅ Covered |
| Poly-T negative (3 T's) | S-009 | ✅ Covered |
| Self-complementarity > 0.3 penalty | C-001 (8bp palindrome, Score=25.625) | ✅ Covered |
| Self-complementarity ≤ 0.3 no penalty | C-001 control (8bp, Score=100) | ✅ Covered |
| Seed GC penalty (low, < 30%) | S-010 (Score=95) | ✅ Covered |
| Seed GC penalty (high, > 80%) | S-011 (Score=95) | ✅ Covered |
| Restriction site penalty | S-001 (Score=95) | ✅ Covered |
| Score clamping to 0 | C-002 (Score=0) | ✅ Covered |
| FullGuideRna scaffold (76nt) | M-006 (exact scaffold match, length=96) | ✅ Covered |
| DesignGuideRnas PAM-based | S-003 (Count=1, verified props) | ✅ Covered |
| DesignGuideRnas multi-PAM | C-005 (Count=3) | ✅ Covered |
| DesignGuideRnas no-PAM | C-004 (empty) | ✅ Covered |
| DesignGuideRnas MinScore filter | C-010 (filtered when MinScore=101) | ✅ Covered |
| SaCas9 system support | C-006 (System.Name, GuideLength) | ✅ Covered |
| Null/empty argument guards | M-005, M-007, M-008, M-009, C-003 | ✅ Covered |
| GuideRnaParameters defaults | S-004 (5 fields) | ✅ Covered |
| GuideRnaParameters custom | S-005 (5 fields) | ✅ Covered |
| Property: GC in [0,100] | Property test, non-vacuous | ✅ Covered |
| Property: Score in [0,100] | Property test, non-vacuous | ✅ Covered |

**No ❌ Missing, no ⚠ Weak, no 🔁 Duplicate code paths remain.**

---

## Sign-Off

- **Author**: QA Architect
- **Date**: 2026-03-03
- **Version**: 3.0
