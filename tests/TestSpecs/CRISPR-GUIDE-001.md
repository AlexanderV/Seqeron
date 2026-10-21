# TestSpec: CRISPR-GUIDE-001 - Guide RNA Design

## Test Unit Identification
- **Test Unit ID**: CRISPR-GUIDE-001
- **Algorithm Group**: MolTools
- **Algorithm Name**: Guide_RNA_Design
- **Canonical Methods**:
  - `CrisprDesigner.DesignGuideRnas(DnaSequence, int, int, CrisprSystemType, GuideRnaParameters?)`
  - `CrisprDesigner.EvaluateGuideRna(string, CrisprSystemType, GuideRnaParameters?)`
  - `CrisprDesigner.GetCutSite(PamSite)` (B08 F3)
  - `CrisprDesigner.CountSelfComplementaryStems(string, int, IEnumerable<string>?)` + `StandardSgRnaBackboneRegion` (B08 F6)
  - `CrisprDesigner.GetGrafMotif(string)` (B08 F7)
- **Result/parameter contracts added by the review**: `GuideRnaCandidate.Context30Mer` / `OnTargetScore` / `SelfComplementaryStems` / `GrafMotif`, `GuideRnaParameters.Ranking` (`GuideRnaRanking`) / `MaxSelfComplementaryStems` (B08 F5–F7)
- **Related Documentation**: [Guide_RNA_Design.md](../../docs/algorithms/MolTools/Guide_RNA_Design.md)

## Evidence Summary

### Primary Sources
| Source | Key Information | URL |
|--------|-----------------|-----|
| Addgene CRISPR Guide | gRNA structure, 20bp spacer, seed sequence (8-10bp at 3'), PAM requirement, Pol III TTTT terminator | https://www.addgene.org/guides/crispr/ |
| Wikipedia: Guide RNA | GC >50% optimal, length 17-24bp standard 20bp | https://en.wikipedia.org/wiki/Guide_RNA |
| Wikipedia: PAM | SpCas9 PAM is 5'-NGG-3', PAM required for cleavage | https://en.wikipedia.org/wiki/Protospacer_adjacent_motif |
| CRISPOR `crispor.py` (maximilianh/crisporWebsite, master) | Cas9 cut "-3bp 5' of the PAM site" (`startFt = start - 3` / `ftSeq + "---"`), Cpf1 cut "after the 18th base on the non-targeted strand" (Zetsche 2015 Fig. 3); `mergeGuideInfo` `guideData.sort(reverse=True, …)` ranking (F3, F5) | https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py |
| CRISPOR `crisporEffScores.py` (same repo) | Rule Set 2 window `trimSeqs(seqs, -24, 6)` (4 nt 5' flank + 20-nt protospacer + PAM + 3 nt 3' flank), `calcAziScore` "can't do Ns" → −1; `getGrafType` TT-/GCC-motifs, NGG only (F5, F7) | https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crisporEffScores.py |
| CHOPCHOP `chopchop.py` (Carl-labhub/chopchop-clone mirror) | `STEM_LEN = 4`, `selfComp` stem count (GC ≥ 0.5, PAM excluded), backbone passed "same strand as guide" and reverse-complemented internally, `-filterSelfCompMax` default −1 (F6); Labun et al. 2016 NAR 44:W272, 2019 NAR 47:W171 | https://github.com/Carl-labhub/chopchop-clone |
| Doench et al. 2016, Nat Biotechnol 34:184 (Rule Set 2 / Azimuth) | On-target score on designed guides; oracle row `nopos_oracle.csv` (F5) | see CRISPR-OFF-001 / Rule Set 2 evidence |
| Graf et al. 2019, Cell Reports 26:1098–1103 | TT- and GCC-motif inefficiency (F7) | — |
| Zetsche et al. 2015, Cell 163:759 | Cas12a staggered cut (F3) | — |

> **Provenance note (WP-SRC, 2026-10-14).** All three sources above are on hosts denied by the
> organization's egress policy (`www.addgene.org`, `en.wikipedia.org` — agent proxy answers 403 to
> CONNECT). **No unit of this campaign has opened them** (superseded for the Addgene row — see the
> update below); they are **unverifiable in this
> environment** and must not be read as checked. They are retained, not deleted: the PAM definitions
> were separately re-confirmed against CRISPOR's own tables (see
> [`CRISPR-PAM-001-Evidence.md`](../../docs/Evidence/CRISPR-PAM-001-Evidence.md)), and the 20 bp
> spacer default is pinned by the executed Doench/Azimuth and CRISPOR references. See
> [`docs/Evidence/B08-Source-Provenance.md`](../../docs/Evidence/B08-Source-Provenance.md) §1, §7c.
>
> **Update (B08 WP3, 2026-10-14) — the Addgene CRISPR Guide is now VERIFIED.** `www.addgene.org`
> answers 200; https://www.addgene.org/guides/crispr/ was fetched 2026-10-14 and, in "Basic Mechanisms
> of CRISPR", reads: *"the seed sequence (8–10 bases at the 3′ end of the gRNA targeting sequence)"* —
> the 8-10 nt range behind `CrisprDesigner.cs` (seed comment in `EvaluateGuideRna`) is confirmed as
> cited. The same page also confirms the ~20-nt spacer + scaffold structure, the SpCas9 `NGG`, SaCas9
> `NNGRRT` and Cas12a `TTTV` PAMs. It does **not** mention a Pol III `TTTT` terminator: that row's
> Addgene attribution remains unverified. The two Wikipedia rows remain unverifiable here.
> Details: [`B08-Source-Provenance.md`](../../docs/Evidence/B08-Source-Provenance.md) §8.

### Evidence-Backed Parameters
| Parameter | Evidence Value | Source | Implementation Value | Status |
|-----------|----------------|--------|---------------------|--------|
| Guide length | 20bp standard (17-24bp range) | Wikipedia Guide RNA | 20bp | ✅ Aligned |
| Optimal GC | >50% | Wikipedia Guide RNA | 40-70% range | ✅ Broader range acceptable per common tools |
| Seed region | 8-10bp at 3' end — verbatim "the seed sequence (8–10 bases at the 3′ end of the gRNA targeting sequence)" (fetched 2026-10-14) | Addgene CRISPR Guide | 10bp (upper bound) | ✅ Aligned (verified) |
| Poly-T termination | TTTT (Pol III terminator) | Addgene CRISPR Guide — **not on the page** (fetched 2026-10-14; unverified) | TTTT detection | ✅ Aligned |
| SpCas9 PAM | 5'-NGG-3' | Addgene, Wikipedia PAM | NGG | ✅ Aligned |
| SaCas9 PAM | NNGRRT | Addgene PAM table | NNGRRT | ✅ Aligned |
| Cas12a PAM | 5'-TTTV-3' | Addgene PAM table | TTTV | ✅ Aligned |
| Scaffold | tracrRNA-derived scaffold | Addgene: sgRNA = spacer + scaffold | SpCas9 76nt scaffold | ✅ Aligned |

### Scoring Model (verified via exact-value tests)
| Penalty | Formula | Verified By |
|---------|---------|-------------|
| GC outside [MinGc, MaxGc] | −(delta × 2) per % | M-001 (0 penalty), M-002 (−80), M-003 (−60) |
| Poly-T (TTTT+) | −20 flat, only when `AvoidPolyT` (F8) | M-004 (score 80), S-008 (score 80), S-017 (flag off → 100) |
| CHOPCHOP self-complementary 4-bp stems | −1 per stem (`SCORE['FOLDING'] = 1`), only when `CheckSelfComplementarity` (F8, E3; replaced the unsourced `selfComp × 30` above 0.3) | C-001 (11 stems → −11), S-018 (flag off → no deduction), C-023 |
| Seed GC outside [30%, 80%] | −5 flat | S-010 (−5), S-011 (−5), M-002 (combined), M-003 (combined) |
| Restriction site (common) | −5 flat | S-001 (score 95) |
| Score floor | max(0, score) | C-002 (clamped to 0) |
| CHOPCHOP 4-bp stems > `MaxSelfComplementaryStems` | 0 (issue in `EvaluateGuideRna`, filter in `DesignGuideRnas`; null = no filter, CHOPCHOP default −1) (F6) | C-022 |
| C-023 | `EvaluateGuideRna_StemPenalty_IsOnePointPerChopchopStem` (5 cases) | ACGT×5, GCGC×10, `ACGGACTAGCCTACGTACGT`, ACGT×3, GCGCGCGC ± `StandardSgRnaBackboneRegion` | stems 10/11/0/2/0 and 10/11/8/2/0 with backbone (chopchop.py `calcSelfComplementarity` run verbatim); Score = unpenalised − stems (E3) | ✅ |
| C-024 | `EvaluateGuideRna_LegacySelfComplementarityScore_DoesNotEnterTheScore` | "GCGCGCGC" | legacy 0.3125, 0 stems, Score=95 (E3) | ✅ |
| C-025 | `RawStringGuideEntryPoints_RejectNonNucleotides` (5 cases) | `X`×20, `N`, `-`, space, `ſ` (U+017F) | `ArgumentException` (ParamName `guideSequence`) from `EvaluateGuideRna`, `CountSelfComplementaryStems` and the backbone list (E3) | ✅ |
| C-026 | `RawStringGuideEntryPoints_FoldCaseAndReadUAsT` | lower-case RNA spacer, RNA scaffold | identical to the DNA upper-case result (E3) | ✅ |
| C-027 | `CalculateMitHitScore_LongS_IsNotFoldedOntoS` | `…GS` vs `…Gſ` | one mismatch at position 19, equal to the `…GA` mismatch score (was 100) (D20/E3) | ✅ |
| Graf 2019 TT-/GCC-motif | 0 (issue only, NGG systems only — CRISPOR warns) (F7) | C-017 |

---

## Test Categories

### MUST Tests (Critical Functionality) — 9 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| M-001 | `EvaluateGuideRna_OptimalGuide_HighScore` | "ACGGACTAGCCTACGTACGT" (55% GC, 0 stems) | Score=100, GC=55, SeedGC=50, Issues=[], Stems=0 (E3: the former ACGT×5 guide carries 10 stems) | ✅ |
| M-002 | `EvaluateGuideRna_LowGcContent_LowerScore` | "AAAAAAAAAAAAAAAAAAAA" (0% GC) | Score=15, GC=0, SeedGC=0, Issues=2 [Low GC, Suboptimal seed GC] | ✅ |
| M-003 | `EvaluateGuideRna_HighGcContent_LowerScore` | "GCGCGCGCGCGCGCGCGCGC" (100% GC) | Score=35, GC=100, SeedGC=100, Issues=2 [High GC, Suboptimal seed GC] | ✅ |
| M-004 | `EvaluateGuideRna_HasPolyT_Penalized` | "ACGTACGTTTTTACGTACGT" (TTTT) | Score=80, GC=40, HasPolyT=true, Issues=1 [TTTT] | ✅ |
| M-005 | `EvaluateGuideRna_EmptyGuide_ThrowsException` | "" | Throws ArgumentNullException | ✅ |
| M-006 | `EvaluateGuideRna_FullGuideRna_IncludesScaffold` | any 20bp guide | FullGuideRna = guide + scaffold (96 chars total) | ✅ |
| M-007 | `DesignGuideRnas_NullSequence_ThrowsException` | null | Throws ArgumentNullException | ✅ |
| M-008 | `DesignGuideRnas_InvalidRegionStart_ThrowsException` | regionStart = −1 | Throws ArgumentOutOfRangeException | ✅ |
| M-009 | `DesignGuideRnas_InvalidRegionEnd_ThrowsException` | regionEnd > len | Throws ArgumentOutOfRangeException | ✅ |

---

### SHOULD Tests (Important Quality) — 19 tests

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

#### Added by the 2026-09 review (B08 F3–F5, F8) — 8 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| S-012 | `DesignGuideRnas_ReverseStrandGuide_CutSiteInRegion_IsDesigned` | "CCA"+ACGT×10, whole region | reverse site Pos=0, TargetStart=3, `GetCutSite`=5; 1 guide, Pos=3, reverse (F3, CRISPOR −3 bp 5' of PAM) | ✅ |
| S-013 | `DesignGuideRnas_ReverseStrandGuide_CutSiteOutsideRegion_IsExcluded` | same, regions [6,30], [0,4], [5,5] | empty, empty, exactly 1 (F3 region boundary) | ✅ |
| S-014 | `GetCutSite_FollowsCrisporConventionOnBothStrands` | SpCas9 ±, Cas12a ± PamSites at 40 | 37 / 45 / 62 / 21 (F3: Cas9 `Position−3` / `Position+pamLen+2`; Cas12a `Position+pamLen+18` / `Position−19`) | ✅ |
| S-015 | `DesignGuideRnas_Output_IsRankedByScoreDescending` | 3 PAMs scoring 100/80/0 | Scores [100, 80, 0] (F5, CRISPOR `mergeGuideInfo` sort) | ✅ |
| S-016 | `DesignGuideRnas_KeepsTheRequestedSystemMetadata` | LbCas12a target | System LbCas12a, GuideLength 24, PamAfterTarget=false, SeedGC=60 (5' seed) (F4) | ✅ |
| S-017 | `EvaluateGuideRna_AvoidPolyTFalse_DisablesThePolyTPenalty` | "ACGTACGTTTTTACGTACGT" | on: 80; off: 100, HasPolyT=true, Issues=[] (F8) | ✅ |
| S-018 | `EvaluateGuideRna_CheckSelfComplementarityFalse_DisablesThePenalty` | GCGC×10 (11 stems) | on: 100−11−5; off: 95, Stems=11, legacy SelfComp=0.275, 1 seed-GC issue (F8, E3) | ✅ |
| S-019 | `DesignGuideRnas_EqualRankingKeys_AreOrderedByPositionDeterministically` (2 cases: `Score`, `OnTargetRuleSet2`) | "CCA"+ACGT×10+"AGG" | 2 tied guides ordered (3, reverse), (23, forward) — deterministic position/strand tie-break (F5) | ✅ |

---

### COULD Tests (Edge Cases) — 22 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| C-001 | `EvaluateGuideRna_SelfComplementary_PenaltyTriggered` | GCGC×10 (11 stems) | Score=100−60−11−5=24, Issues ∋ "11 self-complementary 4-bp stems"; control "GCGCGCGC" (0 stems in the reference, legacy fraction 0.3125) Score=35, no stem issue (E3) | ✅ |
| C-002 | `EvaluateGuideRna_AllT_VeryLowScoreWithMultipleIssues` | "TTTTTTTTTTTTTTTTTTTT" | Score=0 (clamped), GC=0, SeedGC=0, HasPolyT=true, Issues=3 | ✅ |
| C-003 | `EvaluateGuideRna_NullGuide_ThrowsException` | null | Throws ArgumentNullException | ✅ |
| C-004 | `DesignGuideRnas_NoPamInRegion_ReturnsEmpty` | all-A sequence | Returns empty | ✅ |
| C-005 | `DesignGuideRnas_MultiplePams_ReturnsMultipleGuides` | 3 × (20bp+NGG) | Count=3, all Score=100 | ✅ |
| C-006 | `EvaluateGuideRna_SaCas9SystemType_ValidEvaluation` | "ACGTACGTACGTACGTACGT" SaCas9 | Score=100, System.Name="SaCas9", GuideLength=21 | ✅ |
| C-007 | `EvaluateGuideRna_BelowBoundaryGc_HasLowGcIssue` | "AAAAAAAAAAAAGCGCGCAT" (30%) | GC=30, Score=80, Issues=1 [Low GC] | ✅ |
| C-008 | `DesignGuideRnas_EntireSequenceAsRegion_Works` | 27bp seq with PAM | Count=1, Pos=4, Score=100 | ✅ |
| C-009 | `EvaluateGuideRna_AboveBoundaryGc_HasHighGcIssue` | "GCGCGCGCGCGCGCGCAAAA" (80%) | GC=80, Score=80, Issues=1 [High GC] | ✅ |
| C-010 | `DesignGuideRnas_MinScoreFiltering_ExcludesLowScoreGuides` | seq with PAM | MinScore=50 → Count=1; MinScore=101 → empty | ✅ |

#### Added by the 2026-09 review (B08 F5–F7) — 12 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| C-011 | `DesignGuideRnas_ForwardGuide_CarriesRuleSet2ContextAndScore` | Azimuth oracle 30-mer `AAAAAAAAAAAAAAAGCTAACAGCAGGAGT` | `Context30Mer` = the 30-mer, `OnTargetScore` = 0.529037 ± 1e-5 (F5, `nopos_oracle.csv`) | ✅ |
| C-012 | `DesignGuideRnas_ReverseGuide_CarriesTheSameRuleSet2ContextAndScore` | revcomp of the oracle 30-mer | reverse guide Pos=6, same context and score (F5) | ✅ |
| C-013 | `DesignGuideRnas_WithoutA30MerContext_LeavesTheOnTargetScoreNull` | standalone, 3'/5' flank overrun, Cas12a | `Context30Mer` and `OnTargetScore` null (F5) | ✅ |
| C-014 | `DesignGuideRnas_RankingByRuleSet2_OrdersByTheAzimuthScore` | 1 unscorable + 2 scorable loci | scored guides descending by `OnTargetScore`, unscored last (F5, `GuideRnaRanking.OnTargetRuleSet2`) | ✅ |
| C-015 | `GetGrafMotif_MatchesTheCrisporReference` (20 cases) | each `getGrafType` branch | None / TtMotif / GccMotif per CRISPOR (F7) | ✅ |
| C-016 | `GetGrafMotif_ExhaustiveSuffixCensus_MatchesTheReferenceCounts` | all 4096 3'-terminal 6-mers | TT 496, GCC 64, None 3536 (F7) | ✅ |
| C-017 | `EvaluateGuideRna_GrafMotif_IsReportedForNggSystemsWithoutAPenalty` | "ACGTACGTACGTACGTACTT" SpCas9 vs SaCas9 | SpCas9: TtMotif, Score=100, issue; SaCas9: None, no issue (F7) | ✅ |
| C-018 | `CountSelfComplementaryStems_MatchesTheChopchopReference` (13 cases) | guides ± `StandardSgRnaBackboneRegion` | CHOPCHOP `selfComp` counts, e.g. ACGT×5 → 10/10, `ACGGACTAGCCTACGTACGT` → 0/8 (F6) | ✅ |
| C-019 | `CountSelfComplementaryStems_ExhaustiveNineMerCensus_MatchesTheReference` | all 262 144 9-mers | Σ 704 (no backbone) / 41 632 (standard backbone) (F6) | ✅ |
| C-020 | `CountSelfComplementaryStems_BackboneRegions_AreReverseComplementedLikeChopchop` | `ACGGACTAGCCTACGTACGT` | 0 / 8 (guide-strand backbone) / 1 (pre-revcomped backbone) (F6) | ✅ |
| C-021 | `CountSelfComplementaryStems_InvalidArguments_Throw` | null, "", stemLength 0 | ArgumentNullException ×2, ArgumentOutOfRangeException (F6) | ✅ |
| C-022 | `DesignGuideRnas_MaxSelfComplementaryStems_FiltersLikeChopchop` | ACGT×5+"AGG" | 10 stems; limit 10 keeps (inclusive), 9 drops; `EvaluateGuideRna` limit 2 → issue, Score=100 (F6) | ✅ |

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
| Guide RNA Evaluation + Design (50 test methods, 82 test cases) | CrisprDesigner_GuideRNA_Tests.cs |
| Guide RNA Properties (2 tests) | Properties/CrisprProperties.cs |
| Off-Target Analysis | CrisprDesigner_OffTarget_Tests.cs (CRISPR-OFF-001) |
| PAM Detection | CrisprDesigner_PAM_Tests.cs (CRISPR-PAM-001) |

---

## Coverage Summary

| Category | MUST | SHOULD | COULD | Total |
|----------|------|--------|-------|-------|
| Required | 9 | 19 | 22 | 50 |
| Implemented | 9 | 19 | 22 | 50 |
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
| Stem penalty (1 point per CHOPCHOP stem) | C-001, C-023 (oracle-locked stems 10/11/0/2/0, with and without backbone) | ✅ Covered |
| Legacy fraction reported, not scored | C-001 control, C-024 | ✅ Covered |
| Raw-string alphabet (A/C/G/T, U→T, ASCII fold; E3) | C-025, C-026, C-027 | ✅ Covered |
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
| Cut site, both strands / both PAM geometries (F3) | S-012, S-013, S-014 | ✅ Covered |
| System metadata + PAM-before seed orientation kept (F4) | S-016 | ✅ Covered |
| Rule Set 2 context + score, both strands, null cases (F5) | C-011, C-012, C-013 | ✅ Covered |
| Ranking by Score / OnTargetRuleSet2 + tie-break (F5) | S-015, C-014, S-019 | ✅ Covered |
| Rule Set 2 window containing a non-ACGT base → null (F5, CRISPOR "can't do Ns") | — | ⛔ Unreachable via the public API: `DesignGuideRnas` takes a `DnaSequence`, which admits only A/C/G/T, so the guard in `TryGetRuleSet2Context` is defensive only |
| CHOPCHOP stem count, backbone orientation, filter/issue, validation (F6) | C-018 – C-022 | ✅ Covered |
| Graf TT-/GCC-motifs, NGG-only, no deduction (F7) | C-015, C-016, C-017 | ✅ Covered |
| `AvoidPolyT` / `CheckSelfComplementarity` gating (F8) | S-017, S-018 | ✅ Covered |

**No ❌ Missing, no ⚠ Weak, no 🔁 Duplicate code paths remain.**

---

## Sign-Off

- **Author**: QA Architect
- **Date**: 2026-03-03
- **Version**: 3.0
- **Update (B08 WP6/D16, 2026-10-15)**: F3–F8 contracts traced to their tests and sources; S-019 added (tie-break). See `docs/Validation/review-2026-09/B08.md` F3–F8, F47.
