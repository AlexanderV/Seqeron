# MIRNA-TARGET-001: Target Site Prediction — Test Specification

## Test Unit
- **ID**: MIRNA-TARGET-001
- **Area**: MiRNA
- **Canonical methods**: `FindTargetSites(mRna, miRna, minScore)`, scoring via `CalculateTargetScore` (internal); `ScoreTargetSiteContextPlusPlus(mRna, miRna, site, ContextPlusPlusInputs)` (opt-in TargetScan context++, Agarwal 2015 — miRNA+3'UTR-derivable features incl. 3P_score/Min_dist/Len_3UTR/Off6m, plus caller-supplied SPS/TA/Len_ORF/ORF8m)
- **Supporting**: `AlignMiRnaToTarget`, `CreateTargetSite` (internal)

## Canonical Test File
`MiRnaAnalyzer_TargetPrediction_Tests.cs`

## Reference Data
- hsa-let-7a-5p: `UGAGGUAGUAGGUUGUAUAGUU`, seed `GAGGUAG`, seed RC `CUACCUC`
- hsa-miR-21-5p: `UAGCUUAUCAGACUGAUGUUGA`, seed `AGCUUAU`, seed RC `AUAAGCU`

---

## Must Tests (evidence-backed)

### M-001: FindTargetSites — 8mer site detected and scored highest
**Evidence**: Bartel (2009) — 8mer = seed match positions 2-8 + A opposite position 1
- Construct mRNA with exact 8mer site (`CUACCUCA` for let-7a)
- Assert: site found, Type = Seed8mer, SeedMatchLength = 8, score ≥ 0.9

### M-002: FindTargetSites — 7mer-m8 site detected
**Evidence**: Bartel (2009) — 7mer-m8 = match to positions 2-8
- Construct mRNA with full seed RC (pos 2-8) but no trailing A (`CUACCUCG`)
- Assert: site found, Type = Seed7merM8, SeedMatchLength = 7

### M-003: FindTargetSites — 7mer-A1 site detected
**Evidence**: Bartel (2009) — 7mer-A1 = positions 2-7 + A opposite position 1
- Construct mRNA with 6mer core (pos 2-7 RC = `UACCUC`) + downstream A (`UACCUCA`)
- Assert: site found, Type = Seed7merA1, SeedMatchLength = 7

### M-004: FindTargetSites — 6mer site detected
**Evidence**: Bartel (2009) — 6mer = match to positions 2-7 only
- Construct mRNA where only 6-mer (pos 2-7 RC) matches, no A context
- Assert: site found, Type = Seed6mer, SeedMatchLength = 6

### M-005: Score monotonicity — 8mer > 7mer-m8 > 7mer-A1 > 6mer
**Evidence**: Grimson (2007) — efficacy hierarchy is well-established
- Create four mRNAs, each with exactly one site type
- Assert: score(8mer) > score(7mer-m8) > score(7mer-A1) > score(6mer)

### M-006: FindTargetSites — empty/null inputs return empty
**Evidence**: Defensive API contract
- Empty mRNA → empty, null miRNA sequence → empty

### M-007: FindTargetSites — no match returns empty
**Evidence**: Trivial correctness
- mRNA with no seed RC occurrence → empty result

### M-008: FindTargetSites — multiple sites found independently
**Evidence**: Bartel (2009) — each seed match site functions independently
- mRNA with 3 copies of seed RC separated by spacers
- Assert: Count == 3 sites, all classified as 7mer-m8

### M-009: Score range is [0.0, 1.0]
**Evidence**: Implementation contract; Grimson (2007) context scores are bounded
- For all 5 site types (8mer, 7mer-m8, 7mer-A1, 6mer, offset 6mer), assert 0 ≤ score ≤ 1

### M-010: FindTargetSites — minScore filtering works
**Evidence**: API contract
- Use high minScore threshold (0.99); verify low-scoring sites excluded

### M-011: AlignMiRnaToTarget — perfect complementary yields all matches
**Evidence**: Watson-Crick pairing definition
- miRNA = "AAAA", target = "UUUU" → 4 matches, 0 mismatches

### M-012: AlignMiRnaToTarget — G:U wobble detected
**Evidence**: Wobble pair definition (Crick 1966)
- GGGG vs UUUU → 4 wobbles

### M-013: AlignMiRnaToTarget — mismatches counted
**Evidence**: Trivial correctness
- AAAA vs AAAA → 4 mismatches (same bases don't pair)

### M-014: AlignMiRnaToTarget — empty input returns empty duplex
**Evidence**: Defensive API contract

### M-015: Free energy is negative for well-paired duplex
**Evidence**: Thermodynamic principle; Bartel (2009)
- Perfect complement pairing → negative free energy

### M-016: Target site includes alignment string
**Evidence**: API contract
- Found site has non-empty Alignment field

### M-017: DNA input (T) handled equivalently to RNA (U)
**Evidence**: Implementation design — T→U conversion
- mRNA with T bases, same result as U

---

## Should Tests

### S-001: Offset 6mer detected
- mRNA with match to positions 3-8 RC (`CUACCU` for let-7a) → Type = Offset6mer

### S-002: TargetSite record fields populated correctly
- All fields (Start, End, TargetSequence, MiRnaName, Type, Score, FreeEnergy, Alignment) populated

### S-003: Real miRNA (let-7a, miR-21) against constructed targets
- let-7a against 3'UTR-like mRNA → exactly 1 site, 8mer
- miR-21 against constructed 8mer target → exactly 1 site, 8mer

---

### CTX-001..CTX-011: TargetScan context++ scoring (Agarwal et al. 2015) — opt-in `ScoreTargetSiteContextPlusPlus`

**Source (retrieved verbatim this session):** `Agarwal_2015_parameters.txt` (fitted coefficients) and `targetscan_70_context_scores.pl` (feature computation/scaling: `getAgarwalContribution`, `getLocalAU_contribution`, `get_sRNA1_8_contributions`, `getSite8_contribution`, `get3primePairingContribution` + `extractSubseqForAlignment`/`modifySubseqForAlignment`, `getMinDist_weighted_contribution`, `get_len3UTR_weighted_contribution`, `getOffset6merSites`/`getOffset6mer_weighted_contribution`), TargetScan distribution; peer-reviewed model Agarwal V et al. (2015) eLife 4:e05005, doi:10.7554/eLife.05005. Each expected value is hand-derived from the retrieved coefficients × feature values (3P_score raw values cross-checked against the reference perl); a wrong coefficient/formula fails the test.

- **CTX-001 (MUST):** 8mer, let-7a, all-G flanks (len 18) → `ContextScorePartial == -0.7561913315126536` with per-feature: Intercept=-0.589, Local_AU=+0.154608695652174, sRNA8G(8mer)=+0.015, 3P_score=+0.016 (raw 0), Min_dist=-0.049759446106213065, Len_3UTR=-0.2830405810586145, Off6m=-0.020 (one 'CUACCU'). Exact, `Within(1e-9)`.
- **CTX-002 (MUST):** 7mer-m8, let-7a (len 16) → `ContextScorePartial == -0.3087858306402126`; 3P_score=+0.022, Min_dist=-0.031015975380457392, Len_3UTR=-0.15385698397262643, Off6m=-0.011.
- **CTX-003 (MUST):** 7mer-A1, let-7a, Site8 'G' (len 16) → `ContextScorePartial == -0.2673816755691017`; 3P_score=+0.024, Min_dist=-0.022124733327545488, Len_3UTR=-0.12813929518273268, Off6m=0.
- **CTX-004 (MUST):** 6mer, miR-21, mixed flanks (len 16), Site8 'C' → `ContextScorePartial == -0.12128221179295548`; 3P_score=+0.0096, Min_dist=-0.017194033053347654, Len_3UTR=-0.04447703767941017, Off6m=0.
- **CTX-005 (MUST):** sRNA1 non-U branch — synthetic miRNA nt1=G, nt8=G, 8mer → `SRna1Contribution == 0.060` and `SRna8Contribution == 0.015`.
- **CTX-006 (MUST):** non-seed-match site type (Offset6mer) → throws `ArgumentException`.
- **CTX-007 (MUST):** no optional inputs → `OmittedFeatures` contains `SA, PCT, SPS, TA_3UTR, Len_ORF, ORF8m` and does NOT contain `3P_score, Min_dist, Len_3UTR, Off6m` (now computed).
- **CTX-008 (MUST):** 8mer with real 3' supplementary pairing (raw 3P=6, perl-verified) → `ThreePrimePairingContribution == -0.080` (= -0.040×((6-1)/2.5)).
- **CTX-009 (MUST):** 3P raw scores across site types reproduce the perl reference: 7mer-m8 raw 4.5 → contribution -0.040×((4.5-1)/2.5); 7mer-A1 raw 6 → -0.060×((6-1)/2.5).
- **CTX-010 (MUST):** UTR with two 'CUACCU' offset-6mers → `Off6mContribution == -0.040` (raw count 2 × -0.020).
- **CTX-011 (MUST):** caller-supplied `ContextPlusPlusInputs(Sps=-8.0, Ta=3.5, OrfLength=1000, Orf8mCount=2)` on 8mer → SPS=+0.11716577540106952, TA=+0.11424734042553189, Len_ORF=+0.04503626943005184, ORF8m=-0.236; those four drop from `OmittedFeatures`, only `SA, PCT` remain.

---

## Could Tests

### C-001: Context (AU-rich) may influence scoring indirectly via supplementary pairing bonuses

---

## Coverage Classification

**Date**: 2026-03-16
**Test file**: `MiRnaAnalyzer_TargetPrediction_Tests.cs`
**Tests**: 27 (17 Must + 4 Should + 6 parametric variants)
**Result**: 0 missing, 0 weak, 0 duplicate

| Spec | Test Method | Status | Notes |
|------|-----------|--------|-------|
| M-001 | FindTargetSites_8merSite_DetectedAndScoredHighest | ✅ Covered | Exact type + seedMatchLength + score threshold |
| M-002 | FindTargetSites_7merM8Site_Detected | ✅ Covered | Exact type + seedMatchLength |
| M-003 | FindTargetSites_7merA1Site_Detected | ✅ Covered | Exact type + seedMatchLength |
| M-004 | FindTargetSites_6merSite_Detected | ✅ Covered | Exact type + seedMatchLength |
| M-005 | FindTargetSites_ScoreMonotonicity_8mer_GT_7merM8_GT_7merA1_GT_6mer | ✅ Covered | Strict inequality chain |
| M-006 | FindTargetSites_EmptyMrna_ReturnsEmpty + EmptyMiRnaSequence | ✅ Covered | Two tests |
| M-007 | FindTargetSites_NoSeedMatch_ReturnsEmpty | ✅ Covered | |
| M-008 | FindTargetSites_MultipleSites_AllFound | ✅ Covered | Exact count == 3, all 7mer-m8 |
| M-009 | FindTargetSites_AllSites_ScoreInRange (×5) | ✅ Covered | All 5 site types: 8mer, 7mer-m8, 7mer-A1, 6mer, offset 6mer |
| M-010 | FindTargetSites_HighMinScore_FiltersLowScoringSites | ✅ Covered | |
| M-011 | AlignMiRnaToTarget_PerfectComplement_AllMatches | ✅ Covered | Exact counts: 4 matches, 0 mismatches, 0 wobbles |
| M-012 | AlignMiRnaToTarget_GUWobblePairs_Detected | ✅ Covered | Exact counts: 4 wobbles, 0 matches |
| M-013 | AlignMiRnaToTarget_SameBases_AllMismatches | ✅ Covered | Exact count: 4 mismatches |
| M-014 | AlignMiRnaToTarget_EmptyMiRna + EmptyTarget | ✅ Covered | Two tests |
| M-015 | AlignMiRnaToTarget_WellPairedDuplex_NegativeFreeEnergy | ✅ Covered | Thermodynamic principle: ΔG < 0 |
| M-016 | FindTargetSites_FoundSite_HasNonEmptyAlignment | ✅ Covered | |
| M-017 | FindTargetSites_DnaInput_ConvertedToRnaAndMatched | ✅ Covered | Same type + score within 0.001 |
| S-001 | FindTargetSites_Offset6merSite_Detected | ✅ Covered | Exact type + seedMatchLength |
| S-002 | FindTargetSites_FoundSite_AllFieldsPopulated | ✅ Covered | All 8 fields verified |
| S-003 | FindTargetSites_Let7a_RealTargetSequence + MiR21_8merTarget | ✅ Covered | Exact type + count for both |
| CTX-001 | ScoreTargetSiteContextPlusPlus_8merLet7a_GcFlanks_MatchesHandDerivedScore | ✅ Covered | Exact partial CS -0.7561913315126536, Within(1e-9); full per-feature breakdown (incl. 3P/Min_dist/Len_3UTR/Off6m) asserted |
| CTX-002 | ScoreTargetSiteContextPlusPlus_7merM8Let7a_MatchesHandDerivedScore | ✅ Covered | Exact partial CS -0.3087858306402126 + new features |
| CTX-003 | ScoreTargetSiteContextPlusPlus_7merA1Let7a_Site8G_MatchesHandDerivedScore | ✅ Covered | Exact partial CS -0.2673816755691017; Site8 path + new features |
| CTX-004 | ScoreTargetSiteContextPlusPlus_6merMiR21_MixedFlanks_MatchesHandDerivedScore | ✅ Covered | Exact partial CS -0.12128221179295548; non-zero local-AU + new features |
| CTX-005 | ScoreTargetSiteContextPlusPlus_SRna1G_8mer_AddsSRna1GCoefficient | ✅ Covered | sRNA1 non-U branch: sRNA1G=0.060, sRNA8G=0.015 |
| CTX-006 | ScoreTargetSiteContextPlusPlus_NonSeedSiteType_Throws | ✅ Covered | ArgumentException for Offset6mer |
| CTX-007 | ScoreTargetSiteContextPlusPlus_NoOptionalInputs_ReportsResidualFeatures | ✅ Covered | Residual = SA/PCT/SPS/TA/Len_ORF/ORF8m; 3P/Min_dist/Len_3UTR/Off6m NOT residual |
| CTX-008 | ScoreTargetSiteContextPlusPlus_3PrimeSupplementaryPairing_8mer_MatchesScaledRawScore | ✅ Covered | 3P raw 6 (perl-verified) → -0.080 |
| CTX-009 | ScoreTargetSiteContextPlusPlus_3PrimeRawScore_PerlReference | ✅ Covered | 3P raw 4.5 (7m8) and 6 (7A1) → scaled contributions |
| CTX-010 | ScoreTargetSiteContextPlusPlus_TwoOffset6mers_CountedRaw | ✅ Covered | Off6m raw count 2 → -0.040 |
| CTX-011 | ScoreTargetSiteContextPlusPlus_SuppliedInputs_8mer_MatchHandDerivedAndDropFromResidual | ✅ Covered | SPS/TA/Len_ORF/ORF8m supplied → exact contributions + drop from residual |

### Resolved Issues (this cycle)

| Issue | Resolution |
|-------|-----------|
| M-008 weak: `Count >= 3` | Strengthened to `Count == 3` + all sites confirmed as 7mer-m8 |
| M-009 weak: only 2 of 5 site types | Added 3 TestCases: 7mer-A1, 6mer, offset 6mer |
| S-003 weak: `Is.Not.Empty` + `Any(type)` | Strengthened to exact `Count == 1` + `Type == Seed8mer` |
| Duplicate: `MiRnaProperties.FindTargetSites_Score_InRange` | Removed (→ M-009) |
| Duplicate: `MiRnaProperties.FindTargetSites_Positions_WithinBounds` | Removed (→ S-002) |

### Audit Trail (previous cycles)

| Test (old MiRnaAnalyzerTests.cs) | Decision | Resolved |
|------|----------|----------|
| FindTargetSites_PerfectSeedMatch_FindsSite | Replaced with M-001 | ✅ |
| FindTargetSites_8mer_HighestScore | Replaced with M-005 | ✅ |
| FindTargetSites_NoMatch_ReturnsEmpty | Migrated as M-007 | ✅ |
| FindTargetSites_MultipleSites_FindsAll | Migrated as M-008 | ✅ |
| FindTargetSites_EmptySequences_ReturnsEmpty | Migrated as M-006 | ✅ |
| FindTargetSites_IncludesAlignment | Replaced with M-016 | ✅ |
| AlignMiRnaToTarget_PerfectMatch_AllMatches | Migrated as M-011 | ✅ |
| AlignMiRnaToTarget_AllMismatches_NoMatches | Migrated as M-013 | ✅ |
| AlignMiRnaToTarget_WobblePairs_Detected | Migrated as M-012 | ✅ |
| AlignMiRnaToTarget_EmptySequences_ReturnsEmptyDuplex | Migrated as M-014 | ✅ |
| AlignMiRnaToTarget_CalculatesFreeEnergy | Replaced with M-015 | ✅ |
| FullWorkflow_PredictTargets | Replaced with S-003 | ✅ |
| TargetSite_HasAllFields | Replaced with S-002 | ✅ |

---

## Deviations and Assumptions

### Resolved (no remaining deviations from external sources)

All site type definitions and classification logic verified against:
- **Bartel (2009)**: Seed positions 2-8; 8mer/7mer-m8/7mer-A1/6mer hierarchy; A opposite position 1 is downstream on mRNA
- **Grimson (2007)**: Scoring weights 8mer=0.310, 7mer-m8=0.161, 7mer-A1=0.099
- **TargetScan 8.0**: Centered sites removed; offset 6mer (positions 3-8) included as marginal type
- **Agarwal (2015)**: 6mer and offset 6mer have minimal but detectable efficacy

### Documented Simplifications (not deviations)

These are intentional engineering simplifications that do not deviate from the biological definitions in external sources:

| # | Area | Simplification | External Source Reference | Impact |
|---|------|---------------|--------------------------|--------|
| 1 | Energy calculation | Simplified stacking/wobble/mismatch model vs Turner nearest-neighbor thermodynamic parameters | Turner & Mathews (2010) | Approximate ΔG values; does not affect site type classification |
| 2 | Alignment | Positional (parallel index) alignment vs true antiparallel structural alignment | Bartel (2009) | Duplex display is schematic; does not affect seed match detection |
| 3 | Default context scoring | Default `Score` = minimal context adjustment vs full TargetScan context++ model. An **opt-in** `ScoreTargetSiteContextPlusPlus` provides the source-fitted context++ with every miRNA+3'UTR-derivable feature (Intercept, Local_AU, sRNA1/8, Site8, **3P_score, Min_dist, Len_3UTR, Off6m**) plus the data-blocked features (`SPS, TA_3UTR, Len_ORF, ORF8m`) when supplied via `ContextPlusPlusInputs`, all with verbatim Agarwal_2015_parameters.txt coefficients and `targetscan_70_context_scores.pl` formulas. | Grimson (2007), Agarwal (2015) eLife 4:e05005 | Default `Score` unchanged; opt-in CS is a **partial** context++ — `SA` (RNAplfold partition-function accessibility, not MFE-approximated) and `PCT` (multi-species conservation) remain honest residuals, plus `SPS/TA_3UTR/Len_ORF/ORF8m` unless supplied; all reported in `OmittedFeatures`. |
| 4 | Score normalization | Base scores normalized to [0,1] from Grimson (2007) absolute weights | Grimson (2007) | Proportional ratios preserved exactly |

### Previously Fixed Issues (audit trail)

| Issue | Root Cause | Fix | Verified Against |
|-------|-----------|-----|-----------------|
| 6mer and offset 6mer patterns were swapped | `seedRC[0:6]` used for 6mer instead of `seedRC[1:7]` | Corrected: 6mer core = seedRC[1:7] (pos 2-7), offset = seedRC[0:6] (pos 3-8) | Bartel (2009), TargetScan 8.0 |
| 7mer-A1 checked preceding A instead of trailing A | `mrna[pos-1] == 'A'` instead of downstream check | Corrected: A opposite position 1 is downstream on mRNA at `mrna[i+6]` | Bartel (2009) |
| 7mer-A1 used wrong seed substring | Used offset pattern (pos 3-8) instead of 6mer core (pos 2-7) | Corrected: uses 6mer core seedRC[1:7] | Bartel (2009) |
| Scoring weights were arbitrary | 1.0/0.9/0.85/0.7/0.6 with no source | Corrected: 1.0/0.52/0.32/0.15/0.10 proportional to Grimson (2007) | Grimson (2007) |
| Guard condition too restrictive | `pos + 8 > mrna.Length` rejected valid 6mer sites | Corrected: scan uses `i <= mrna.Length - 6` | N/A (logic error) |
| Tests fitted to buggy code | mRNA constructions matched wrong patterns | All tests rewritten with correct biological patterns | Bartel (2009) |
