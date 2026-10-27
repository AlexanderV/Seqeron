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

**B08 round 3 G5.** `Score = clamp(100 - penalty, 0, 100)` where `penalty` is CHOPCHOP's own ascending penalty sum (`chopchop.py`, sha256 `3c7a4344ccf2a77e70b638146989d7b9569665216bc15b476177db37630d5960`). Oracle: the `calcGCContent` body, `selfComp`, `SCORE`, `GC_LOW = 40`, `GC_HIGH = 70` copied verbatim into Python 3 + Biopython 1.88 and run on each guide (`EvaluateGuideRna_Score_MatchesChopchopPenaltyOracle`, 16 rows). CHOPCHOP ships the GC term off in every mode (`SCORE_GC = False`; `-c` is `store_false`), so `GuideRnaParameters.ScoreGcContent` (default `true`, a declared deviation) = `false` reproduces CHOPCHOP's default. TTTT / seed GC / cloning site are issue-only (CRISPOR `crispor.py` L3026-3054 prints warnings, no score change). `MinScore` is not from any tool.

| Term | Formula | Verified By |
|---------|---------|-------------|
| GC outside [MinGc, MaxGc] (edges inclusive), when `ScoreGcContent` | +300 flat (`SCORE['CRISPR_BAD_GC']`), i.e. Score 0 | M-001 (in window), M-002 (0 %), M-003 (100 %), S-002 (20 %), C-007 (30 %), C-009 (80 %), S-006/S-007 (edges), oracle rows 20/25/35/40/45/50/60/70/75/80/100 % |
| Poly-T (TTTT+) | issue only, when `AvoidPolyT` (F8; was -20 until G5) | M-004 (95), S-008 (96), S-017 (flag only changes the issue) |
| CHOPCHOP self-complementary 4-bp stems | +1 per stem (`SCORE['FOLDING'] = 1`), only when `CheckSelfComplementarity` (F8, E3) | C-001 (11 stems -> Score 0 with the GC term, 89 without), S-018 (flag off -> no deduction), C-023 |
| Seed GC outside [30%, 80%] | issue only (was -5 until G5) | S-010 (99), S-011 (100) |
| Restriction site (common) | issue only (was -5 until G5) | S-001 (97) |
| Score range | `Math.Clamp(100 - penalty, 0, 100)` | C-002 (GC 0 % -> 0) |
| CHOPCHOP 4-bp stems > `MaxSelfComplementaryStems` | 0 (issue in `EvaluateGuideRna`, filter in `DesignGuideRnas`; null = no filter, CHOPCHOP default −1) (F6) | C-022 |
| C-023 | `EvaluateGuideRna_StemPenalty_IsOnePointPerChopchopStem` (5 cases) | ACGT×5, GCGC×10, `ACGGACTAGCCTACGTACGT`, ACGT×3, GCGCGCGC ± `StandardSgRnaBackboneRegion` | stems 10/11/0/2/0 and 10/11/8/2/0 with backbone (chopchop.py `calcSelfComplementarity` run verbatim); Score = unpenalised − stems (E3) | ✅ |
| C-024 | `EvaluateGuideRna_SelfComplementarityScore_IsTheChopchopStemFraction_AndDoesNotEnterTheScore` | "GCGCGCGC" | `SelfComplementarityScore` = 0/(8−4) = 0 (was the legacy 0.3125), 0 stems, Score=100 (E3, G5, G13) | ✅ |
| C-031 | `EvaluateGuideRna_SelfComplementarityScore_MatchesOracleStemFraction` (6 cases), `..._FollowsBackboneStems` | ACGT×5, GCGC×10, `ACGGACTAGCCTACGTACGT`, ACGT×3, GCGCGCGC, ACGT; standard backbone | stems/(len−4): 10/16, 11/16, 0, 2/8, 0, 0 (no window, never NaN); with backbone 8/16 = 0.5 (G13; stems from chopchop.py `selfComp` run verbatim) | ✅ |
| C-025 | `RawStringGuideEntryPoints_RejectNonNucleotides` (5 cases) | `X`×20, `N`, `-`, space, `ſ` (U+017F) | `ArgumentException` (ParamName `guideSequence`) from `EvaluateGuideRna`, `CountSelfComplementaryStems` and the backbone list (E3) | ✅ |
| C-026 | `RawStringGuideEntryPoints_FoldCaseAndReadUAsT` | lower-case RNA spacer, RNA scaffold | identical to the DNA upper-case result (E3) | ✅ |
| C-027 | `CalculateMitHitScore_LongS_IsRejectedNotFoldedOntoS` | `…GS` vs `…Gſ` | was: one mismatch at position 19 (D20/E3); since G2 `ſ` and `S` are not nucleotides, `ArgumentException` (`offTarget20`), never folded onto `S` | ✅ |
| C-028 | `FindOffTargets_RnaSpellingOfTheGuide_EqualsDnaSpelling_AndDropsTheOnTarget`, `SpecificityAggregates_RnaSpellingOfTheGuide_EqualsDnaSpelling` | RNA / lower-case spelling of a guide | identical sites, PAMs and heuristic / MIT / alt-PAM MIT / CFD scores as the DNA spelling (G2) | ✅ |
| C-029 | `RawStringOffTargetEntryPoints_RejectNonNucleotides` (4), `RawStringOffTargetEntryPoints_KeepTheirGuardOrder` | `N`x20, `X`x20, `-`, `ſ` | `ArgumentException` (`guideSequence`) from all 9 raw-guide entry points; null/empty, genome, range and length guards first (G2) | ✅ |
| C-030 | `GetGrafMotif_RnaSpellingAndCase_EqualDnaSpelling`, `MitHitScore_RejectsNonNucleotides_AndReadsUAsT` | U / lower-case spelling; `N`x20 | same motif / hit score as DNA; `N` rejected (`guide20` / `offTarget20`); CFD / SaCas9 hit scores still reject junk (G2) | ✅ |
| Graf 2019 TT-/GCC-motif | 0 (issue only, NGG systems only — CRISPOR warns) (F7) | C-017 |

---

## Test Categories

### MUST Tests (Critical Functionality) — 9 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| M-001 | `EvaluateGuideRna_OptimalGuide_HighScore` | "ACGGACTAGCCTACGTACGT" (55% GC, 0 stems) | Score=100, GC=55, SeedGC=50, Issues=[], Stems=0 (E3: the former ACGT×5 guide carries 10 stems) | ✅ |
| M-002 | `EvaluateGuideRna_LowGcContent_LowerScore` | "AAAAAAAAAAAAAAAAAAAA" (0% GC) | Score=0 (flat 300 clamped), GC=0, SeedGC=0, Issues=2 [Low GC, Suboptimal seed GC] | ✅ |
| M-003 | `EvaluateGuideRna_HighGcContent_LowerScore` | "GCGCGCGCGCGCGCGCGCGC" (100% GC) | Score=0 (311 clamped), GC=100, SeedGC=100, Issues=3 [High GC, Suboptimal seed GC, 11 stems] | ✅ |
| M-004 | `EvaluateGuideRna_HasPolyT_Penalized` | "ACGTACGTTTTTACGTACGT" (TTTT) | Score=95 (5 stems), GC=40, HasPolyT=true, Issues=2 [TTTT, 5 stems] | ✅ |
| M-005 | `EvaluateGuideRna_EmptyGuide_ThrowsArgumentException` | "" | Throws ArgumentException (null: ArgumentNullException; B08 round 4 H19) | ✅ |
| M-006 | `EvaluateGuideRna_FullGuideRna_IncludesScaffold` | any 20bp guide | FullGuideRna = guide + scaffold (96 chars total) | ✅ |
| M-007 | `DesignGuideRnas_NullSequence_ThrowsException` | null | Throws ArgumentNullException | ✅ |
| M-008 | `DesignGuideRnas_InvalidRegionStart_ThrowsException` | regionStart = −1 | Throws ArgumentOutOfRangeException | ✅ |
| M-009 | `DesignGuideRnas_InvalidRegionEnd_ThrowsException` | regionEnd > len | Throws ArgumentOutOfRangeException | ✅ |

---

### SHOULD Tests (Important Quality) — 19 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| S-001 | `EvaluateGuideRna_RestrictionSite_Penalized` | "ACGTGAATTCACGTACGTAC" (EcoRI) | Score=97 (3 stems), GC=45, Issues=2 [restriction site, 3 stems] | ✅ |
| S-002 | `EvaluateGuideRna_CalculatesSeedGc` | "AAAAAAAAAAAAACGTACGT" | GC=20, SeedGC=40, Score=0 (flat 300), Issues=1 [Low GC] | ✅ |
| S-003 | `DesignGuideRnas_WithPamInRegion_ReturnsGuides` | seq with AGG PAM | Count=1, Pos=24, Score=100, Forward=true | ✅ |
| S-004 | `GuideRnaParameters_Default_HasValidValues` | Default | MinGC=40, MaxGC=70, MinScore=50, AvoidPolyT=true, CheckSelfComp=true | ✅ |
| S-005 | `GuideRnaParameters_CustomValues_Respected` | Custom(30,80,40,false,false) | All 5 fields preserved | ✅ |
| S-006 | `EvaluateGuideRna_BoundaryGc40Percent_NotPenalized` | "AAAAAAAAAAAAGCGCGCGC" (40%) | Score=100, GC=40, SeedGC=80, Issues=[] | ✅ |
| S-007 | `EvaluateGuideRna_BoundaryGc70Percent_NotPenalized` | "GCGCGCGCGCGCGCAAAAAA" (70%) | Score=100, GC=70, SeedGC=40, Issues=[] | ✅ |
| S-008 | `EvaluateGuideRna_ExactlyFourTs_TriggersPolyT` | "ACGTACGTACGATTTTACGT" (4 T's) | Score=96 (4 stems), HasPolyT=true, Issues=2 [TTTT, 4 stems] | ✅ |
| S-009 | `EvaluateGuideRna_ThreeConsecutiveTs_NoPolyT` | "ACGTACGTACGTACGTTTAC" (3 T's) | Score=100, HasPolyT=false, Issues=[] | ✅ |
| S-010 | `EvaluateGuideRna_SeedGcLow_Penalized` | "GCGCGCGCGCAAAAAAAAAA" | GC=50, SeedGC=0, Score=99 (1 stem), Issues=2 [Suboptimal seed GC, 1 stem] | ✅ |
| S-011 | `EvaluateGuideRna_SeedGcHigh_Penalized` | "AAAAAAAAAAGGGGGGGGGG" | GC=50, SeedGC=100, Score=100 (0 stems), Issues=1 [Suboptimal seed GC] | ✅ |

#### Added by the 2026-09 review (B08 F3–F5, F8) — 8 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| S-012 | `DesignGuideRnas_ReverseStrandGuide_CutSiteInRegion_IsDesigned` | "CCA"+ACGT×10, whole region | reverse site Pos=0, TargetStart=3, `GetCutSite`=5; 1 guide, Pos=3, reverse (F3, CRISPOR −3 bp 5' of PAM) | ✅ |
| S-013 | `DesignGuideRnas_ReverseStrandGuide_CutSiteOutsideRegion_IsExcluded` | same, regions [6,30], [0,4], [5,5] | empty, empty, exactly 1 (F3 region boundary) | ✅ |
| S-014 | `GetCutSite_FollowsCrisporConventionOnBothStrands` | SpCas9 ±, Cas12a ± PamSites at 40 | 37 / 45 / 62 / 21 (F3: Cas9 `Position−3` / `Position+pamLen+2`; Cas12a `Position+pamLen+18` / `Position−19`) | ✅ |
| S-014b | `GetCutSite_EveryMember_FollowsItsConvention_OnBothStrands` (`CrisprDesigner_Systems_Tests`, 9 cases) | a PamSite at 40 per system, both strands | PAM-after: 37 / 40+pamLen+2; Cas12a family: 40+pamLen+18 / 21; CasX: 40+4+13=57 / 26 (H38, H51; Liu 2019) | ✅ |
| S-014c | `DesignGuideRnas_EveryMember_IncludesAGuideExactlyWhenItsCutSiteIsInTheRegion` (9 cases) | region = the cut base / ending one before / starting one after | exactly 1 guide / none / none, per system (H51) | ✅ |
| S-014d | `GetSystem_EveryMember_IsPinned`, `GetSystem_CoversEveryEnumValue_AndRejectsAnUnknownOne` | all nine systems | Name, PAM, GuideLength, PamAfterTarget, `CutOffset` (3 / 18 / 13) and `Description` verbatim; `(CrisprSystemType)999` → `ArgumentException` "Unknown CRISPR system" (H44, H50) | ✅ |
| S-015 | `DesignGuideRnas_Output_IsRankedByScoreDescending` | 4 PAMs: 100 (50 % GC, 0 stems) / 95 (5 stems) / 90 (10 stems) / all-T (GC outside window) | Scores [100, 95, 90]; the all-T guide is excluded by the GC filter (F5, CRISPOR `mergeGuideInfo` sort; H8) | ✅ |
| S-016 | `DesignGuideRnas_KeepsTheRequestedSystemMetadata` | LbCas12a target | System LbCas12a, GuideLength 23, PamAfterTarget=false, SeedGC=60 (5' seed) (F4; 24 before H11) | ✅ |
| S-017 | `EvaluateGuideRna_AvoidPolyTFalse_DisablesThePolyTIssue_ScoreUnchanged` | "ACGTACGTTTTTACGTACGT" | on: 95 with a TTTT issue; off: 95, HasPolyT=true, no TTTT issue (F8, G5) | ✅ |
| S-018 | `EvaluateGuideRna_CheckSelfComplementarityFalse_DisablesThePenalty` | GCGC×10 (11 stems) | on: 100−11−5; off: 95, Stems=11, legacy SelfComp=0.275, 1 seed-GC issue (F8, E3) | ✅ |
| S-019 | `DesignGuideRnas_EqualRankingKeys_AreOrderedByPositionDeterministically` (2 cases: `Score`, `OnTargetRuleSet2`) | "CCA"+ACGT×10+"AGG" | 2 tied guides ordered (3, reverse), (23, forward) — deterministic position/strand tie-break (F5) | ✅ |

---

### COULD Tests (Edge Cases) — 22 tests

| ID | Test | Input | Exact Expected | Status |
|----|------|-------|----------------|--------|
| C-001 | `EvaluateGuideRna_SelfComplementary_PenaltyTriggered` | GCGC×10 (11 stems) | Score=0 (311 = 300 GC + 11 stems; 89 with `ScoreGcContent=false`), Issues ∋ "11 self-complementary 4-bp stems"; control "GCGCGCGC" (0 stems in the reference, legacy fraction 0.3125) Score=0 (GC 100 %), no stem issue (E3, G5) | ✅ |
| C-002 | `EvaluateGuideRna_AllT_VeryLowScoreWithMultipleIssues` | "TTTTTTTTTTTTTTTTTTTT" | Score=0 (clamped), GC=0, SeedGC=0, HasPolyT=true, Issues=3 | ✅ |
| C-003 | `EvaluateGuideRna_NullGuide_ThrowsException` | null | Throws ArgumentNullException | ✅ |
| C-004 | `DesignGuideRnas_NoPamInRegion_ReturnsEmpty` | all-A sequence | Returns empty | ✅ |
| C-005 | `DesignGuideRnas_MultiplePams_ReturnsMultipleGuides` | 3 × (20bp+NGG) | Count=3, all Score=100 | ✅ |
| C-006 | `EvaluateGuideRna_SaCas9SystemType_ValidEvaluation` | "ACGTACGTACGTACGTACGT" SaCas9 | Score=100, System.Name="SaCas9", GuideLength=21 | ✅ |
| C-007 | `EvaluateGuideRna_BelowBoundaryGc_HasLowGcIssue` | "AAAAAAAAAAAAGCGCGCAT" (30%) | GC=30, Score=0 (flat 300), Issues=1 [Low GC] | ✅ |
| C-008 | `DesignGuideRnas_EntireSequenceAsRegion_Works` | 27bp seq with PAM | Count=1, Pos=4, Score=100 | ✅ |
| C-009 | `EvaluateGuideRna_AboveBoundaryGc_HasHighGcIssue` | "GCGCGCGCGCGCGCGCAAAA" (80%) | GC=80, Score=0 (300 + 7 stems, clamped), Issues=2 [High GC, 7 stems] | ✅ |
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
| GC penalty (low) | M-002 (Score=0) | ✅ Covered |
| GC penalty (high) | M-003 (Score=0) | ✅ Covered |
| GC boundary 40% (inclusive, no penalty) | S-006 (Score=100) | ✅ Covered |
| GC boundary 70% (inclusive, no penalty) | S-007 (Score=100) | ✅ Covered |
| GC below boundary 30% | C-007 (Score=0) | ✅ Covered |
| GC above boundary 80% | C-009 (Score=0) | ✅ Covered |
| Poly-T detection (4+ T's) | M-004, S-008 | ✅ Covered |
| Poly-T negative (3 T's) | S-009 | ✅ Covered |
| Stem penalty (1 point per CHOPCHOP stem) | C-001, C-023 (oracle-locked stems 10/11/0/2/0, with and without backbone) | ✅ Covered |
| Legacy fraction reported, not scored | C-001 control, C-024 | ✅ Covered |
| Raw-string alphabet (A/C/G/T, U→T, ASCII fold; E3) | C-025, C-026, C-027 | ✅ Covered |
| Seed GC flag (low, < 30%; issue only) | S-010 (Score=99) | ✅ Covered |
| Seed GC flag (high, > 80%; issue only) | S-011 (Score=100) | ✅ Covered |
| Restriction site flag (issue only) | S-001 (Score=97) | ✅ Covered |
| Score clamping to 0 | C-002 (Score=0), M-002, M-003 |
| CHOPCHOP oracle table (GC 20-100 %, TTTT, seed GC, EcoRI, switch on/off) | `EvaluateGuideRna_Score_MatchesChopchopPenaltyOracle` (16 rows), `EvaluateGuideRna_TtttSeedGcCloningSite_AreIssuesOnly`, `GuideRnaParameters_ScoreGcContent_DefaultAndCustomWindow` | ✅ Covered (G5) | ✅ Covered |
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
| Cut site, both strands / both PAM geometries (F3); every one of the nine systems (H38, H51) | S-012, S-013, S-014, S-014b, S-014c | ✅ Covered |
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

## Round 4 contracts and fixes (B08 H8-H12, H14-H23, H41)

Tests in `CrisprDesigner_Contracts_Tests.cs` unless named otherwise. Oracles are named per row.

| ID | Test | Pins | Oracle |
|---|---|---|---|
| H8 | `DesignGuideRnas_GcWindow_FiltersWhetherOrNotTheScoreTermIsOn_Spcas9` (0..20), `..._FiltersFor21NtSaCas9` (0..21), `DesignGuideRnas_FullGcWindow_ListsEveryGuide_AndEvaluateOnlyReports` | the GC window is an independent filter; 0..100 = none | CHOPCHOP chopchop.py L1911-1914: pass sets 8..14 of 20, 9..14 of 21 |
| H9/H41 | `FullGuideRna_FollowsTheGuidesOwnSystem` (9 systems), `SgRnaScaffolds_AreTheOnesCrisporDerivesItsConstantPrimersFrom` | per-system scaffold, `null` for Cas12a/CasX | revcomp of CRISPOR's SpCas9 `tracrPrimer` = scaffold + TTTT; SaCas9 template = scaffold + TTTTTT |
| H10 | `Doench2014_ReferenceExample1/2_MatchesTheExecutedReference`, all-A row | executed `calcDoenchScore` values at 1e-9, docstring literals rejected | `doenchScore.py` executed |
| H11 | `FindPamSites_AllCas12aAliases_ReportTheSameSites`, `FindPamSites_LbCas12a_DetectsTTTV_With23bpTarget` | one Cas12a spacer length | CRISPOR `GUIDELEN = 23` for every Cpf1 PAM |
| H12 | three `FindPamSites_*N*` tests in `CrisprDesigner_PAM_Tests.cs` | CRISPOR `patMatch` / `flankSeqIter(doFilterNs=True)` | genuine `findAllPams` + `flankSeqIter`, 0/0 differential |
| H14 | `CfdSpecificity_List_RefusesNegativeAndNonFiniteScores_WithAMessageThatIsTrue`, `CfdSpecificity_ParkinsonNormalized_ZeroSum_BlamesTheArgumentTheCallerPassed` | negative / NaN sums refused with a true message; the zero-sum error blames `aggregation` for the genome overloads | `calcCfdGuideScore(-2) = -5000.0` is meaningless input |
| H15 | `MitSpecificity_List_RefusesAScoreThatIsNotFiniteAndNonNegative`, `..._KeepsTheDocumentedRangeForValidScores` | the documented [0, 100] | `calcMitGuideScore` is only fed non-negative sums |
| H16 | `GenomeAggregates_FailIdenticallyOnTheSameMalformedCall` (7 malformed calls x 4 aggregates) | one argument prologue | FindOffTargets' guard order + the model refusals |
| H17 | `EvaluateGuideRna_WrongLength_IsRefused_NotScoredAsPerfect`, `..._AcceptsExactlyTheGuideLengthOfEverySystem` | guide length = system length | FindOffTargets' own check |
| H18 | `GuideRnaParameters_Invalid_IsRefusedByDesignAndEvaluation` (7 cases), `..._ParameterlessConstructor_IsTheDocumentedDefault` | validated parameters; `-1` refused | CHOPCHOP `filterSelfCompMax = -1` is "none" |
| H19 | `StringGuideEntryPoints_NullIsAMissingArgument_EmptyIsAWrongValue` (12 entry points) | null -> `ArgumentNullException`, empty -> `ArgumentException` | AzimuthRuleSet2 convention; model_comparison.py L523-L524 |
| H20/H53 | `FindPamSites_UnknownSystem_ThrowsAtTheCall_OnBothOverloads`, `FindPamSites_StringOverload_NullEmptyAndTooShort_YieldNoSites` | eager system validation; null/empty/shorter-than-PAM = no sites | PAM_Site_Detection.md section 3.3 |
| H21 | `RuleSet2_AminoAcidCutPosition_NegativeIsRefused_UnlessThePercentPeptideSentinelIgnoresIt` | documented behaviour of the negative cut position | model_comparison.py L539/L563 |
| H22 | `RuleSet2_NonNggPam_IsScored_AsCrisporDoesWithoutPamAudit` (5), `DesignedGuides_OfTheNagAndNgaSystems_CarryARuleSet2Score` | NAG / NGA contexts scored | tree walker over the original pickles; 300 random contexts (281 non-GG), both models, max diff 0.0 |
| H23 | `ModelBlob_Malformed_FailsAsInvalidOperationException` (8 malformed blobs) | overflow-safe loader, index checks | the documented binary format |
