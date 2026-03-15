# MIRNA-TARGET-001: Target Site Prediction — Test Specification

## Test Unit
- **ID**: MIRNA-TARGET-001
- **Area**: MiRNA
- **Canonical methods**: `FindTargetSites(mRna, miRna, minScore)`, scoring via `CalculateTargetScore` (internal)
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
| 3 | Context scoring | Minimal context adjustment (supplementary pairing bonus, mismatch penalty) vs full TargetScan context++ model | Grimson (2007), Agarwal (2015) | Scores approximate; relative ordering preserved |
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
