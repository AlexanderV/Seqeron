# TestSpec: PRIMER-DESIGN-001

## Primer Pair Design

**Test Unit ID:** PRIMER-DESIGN-001
**Area:** MolTools
**Status:** ☑ Complete
**Last Updated:** 2026-10-01
**Total Tests:** 160 (canonical + smoke + mutation-killing)

---

## Evidence Summary

### Authoritative Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia: Primer (molecular biology)](https://en.wikipedia.org/wiki/Primer_(molecular_biology)) | Encyclopedia | 18-24 bp length, 40-60% GC, Tm 50-60°C, primer pairs within 5°C |
| [Addgene: How to Design a Primer](https://www.addgene.org/protocols/primer-design/) | Protocol Guide | Length 18-24, GC 40-60%, Tm 50-60°C, pairs within 5°C, avoid complementary regions |
| [Primer3 Manual (v2.6.1)](https://primer3.org/manual.html) | Software Documentation | PRIMER_MIN_SIZE=18, PRIMER_MAX_SIZE=27, PRIMER_OPT_SIZE=20, PRIMER_MIN_TM=57, PRIMER_OPT_TM=60, PRIMER_MAX_TM=63, PRIMER_MIN_GC=20, PRIMER_MAX_GC=80, PRIMER_MAX_POLY_X=5, PRIMER_PAIR_MAX_DIFF_TM=100.0 |
| SantaLucia (1998) PNAS 95:1460-65 | Research Paper | Table 2 NN + terminal-initiation parameters — the Primer3-default primer Tm used by EvaluatePrimer/DesignPrimers (the scale of the Primer3 57–63 °C window) |
| Primer3 source `oligotm.c`, `libprimer3.cc` (primer3-org/primer3 main) | Reference implementation | `seqtm`/`oligotm` (SantaLucia Tm + salt correction, `divalent_to_monovalent` 120·√(Mg−dNTP), `long_seq_tm` > 36 nt); `p_obj_fn` penalty; `choose_pair_or_triple` pair search; `primer_rec_comp` / `compare_primer_pair` ordering |
| primer3-py 2.3.1 (`calc_tm`, `design_primers`) | Reference implementation | Numeric oracle for Tm, penalties and selected pairs (values locked in the tests below) |

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
| Pair Tm Difference | ≤ 5°C (unrounded Tm) | Wikipedia, Addgene | Exact match (Primer3 default PRIMER_PAIR_MAX_DIFF_TM=100.0 is unlimited; 5°C is the standard lab guideline) |
| Tm model | Primer3 default (SantaLucia 1998 NN, SantaLucia salt, 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM) | Primer3 PRIMER_TM_FORMULA=1, PRIMER_SALT_CORRECTIONS=1 | Bit-identical to primer3.calc_tm |
| Reaction conditions | `PrimerParameters.MonovalentMillimolar` / `DivalentMillimolar` / `DntpMillimolar` / `DnaConcentrationNanomolar` (null = 50 mM / 1.5 mM / 0.6 mM / 50 nM) for primer Tm, ntthal structure, pair compl and product Tm; internal oligo `Primer3ProbeSettings` conditions | Primer3 PRIMER_SALT_MONOVALENT / _SALT_DIVALENT / PRIMER_DNTP_CONC / PRIMER_DNA_CONC (`p_args`, `create_thal_arg_holder`, `long_seq_tm`), PRIMER_INTERNAL_* (`o_args`) | 2500/2500 random-condition templates identical to primer3-py design_primers (audit round 3, A3-1) |
| GC optimum | `PrimerParameters.OptimalGcPercent`, `Primer3ProbeSettings.OptGcPercent` (+ `WeightGcPercentGt/Lt`); null = 50 % | Primer3 PRIMER_OPT_GC_PERCENT / PRIMER_INTERNAL_OPT_GC_PERCENT (manual default 50; code: undefined, GC weights require it) | Exact match whenever Primer3 accepts the settings |
| Ranking / pair selection | Lowest Primer3 pair objective (`obj_fn`, default Σ per-primer penalty; PRIMER_PAIR_WT_* configurable) over all compatible pairs; Primer3 tie-break; PRIMER_NUM_RETURN ranked pairs (`DesignPrimerPairs`) | Primer3 `choose_pair_or_triple`, `characterize_pair`, `obj_fn`, `compare_primer_pair` | 1000/1000 random templates identical to primer3-py design_primers (ranks 0–4) with Primer3 defaults |
| Product size range | PRIMER_PRODUCT_SIZE_RANGE, default 100–300 bp, ranges tried in order | Primer3 `pr_set_default_global_args_1` (pr_min/pr_max = 100/300), `choose_pair_or_triple` | Exact match (replaces the former ±200 bp flanks) |
| Internal oligo | PRIMER_PICK_INTERNAL_OLIGO: lowest-penalty Primer3 internal oligo strictly between the primers | Primer3 `choose_internal_oligo` | primer3-py PRIMER_INTERNAL_k_* |
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
| `DesignPrimers(template, start, end, params, pairOptions)` | PrimerDesigner | Canonical | O(c log c + p) |
| `DesignPrimerPairs(template, start, end, params, pairOptions)` | PrimerDesigner | Canonical (PRIMER_NUM_RETURN) | O(c log c + p) |
| `CalculateProductMeltingTemperaturePrimer3(product, …)` | PrimerDesigner | Helper (`long_seq_tm`) | O(n) |
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
| M14 | Primer3-default Tm equals primer3.calc_tm (default + non-default conditions, symmetric, > 36 nt) | Tm scale of the Primer3 window | Primer3 oligotm.c, primer3-py |
| M15 | DesignPrimers returns the lowest-penalty compatible pair (primer3-py design_primers), incl. when the individually best primers are Tm-incompatible | Pair selection | Primer3 libprimer3.cc |
| M16 | DesignPrimerPairs with Primer3 defaults returns primer3-py's PRIMER_LEFT/RIGHT_k, PRIMER_PAIR_k_PENALTY / _PRODUCT_TM / _COMPL_ANY_TH / _COMPL_END_TH / _PRODUCT_SIZE for k = 0..4 | PRIMER_NUM_RETURN, product Tm | Primer3 choose_pair_or_triple, long_seq_tm |
| M17 | Default PRIMER_PRODUCT_SIZE_RANGE is 100–300; custom ranges are tried in order; product Tm limits and PRIMER_PAIR_MAX_DIFF_TM are configurable | Search region | Primer3 pr_set_default_global_args_1, characterize_pair |
| M18 | Non-default PRIMER_PAIR_WT_* weights give primer3-py's PRIMER_PAIR_k_PENALTY | Pair objective | Primer3 obj_fn |
| M19 | PRIMER_PICK_INTERNAL_OLIGO picks primer3-py's PRIMER_INTERNAL_k (inside the product, not overlapping the primers; PRIMER_PAIR_WT_IO_PENALTY) | Internal oligo | Primer3 choose_internal_oligo |
| M20 | Primer3 data-control errors throw (weight without optimum, max size > min product, NUM_RETURN < 1, target outside included region) | Validation | Primer3 _pr_data_control |
| M23 | PRIMER_GC_CLAMP / PRIMER_MAX_END_GC / PRIMER_MAX_END_STABILITY reject exactly the primers primer3-py rejects (left/right primers) and give its ranked pairs; illegal values throw; defaults (0 / 5 / 100) inactive; the deprecated `Check3PrimeStability` (−9 kcal/mol gate, unreachable) has no effect and `Avoid3PrimeGC` keeps its library rule | 3′-end checks | Primer3 calc_and_check_oligo_features, _pr_data_control |
| M24 | PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE (+ PRIMER_MIN_THREE_PRIME_DISTANCE) exclude primers whose 3′ ends are too close to those of already selected pairs, giving primer3-py's ranked pairs; −1 (default) reuses primers; < −1 throws | Pair set diversity | Primer3 choose_pair_or_triple, left/right_oligo_in_pair_overlaps_used_oligo |
| M22 | Non-default PRIMER_SALT_MONOVALENT / _DIVALENT / PRIMER_DNTP_CONC / PRIMER_DNA_CONC / PRIMER_OPT_GC_PERCENT (+ PRIMER_INTERNAL_* counterparts) give primer3-py's Tm, structure values, product Tm, penalties and ranking; illegal conditions throw; null = Primer3 defaults | Reaction conditions | Primer3 libprimer3.c p_args/o_args, _pr_data_control |
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
| S4 | ✅ | `EvaluatePrimer_OptimalPrimer_HasHighScore`, `_SuboptimalLength_ScoreVaries`, `EvaluatePrimer_Penalty_IsPrimer3PerPrimerPenalty` | Informational Score; Primer3 penalty locked to primer3-py |
| M14 | ✅ | `CalculateMeltingTemperaturePrimer3_DefaultConditions_MatchesPrimer3CalcTm` (7 cases), `_NonDefaultConditions_…`, `_InvalidInput_NaNOrThrows`, `EvaluatePrimer_NonAcgtBase_TmNotComputableAndInvalid` | primer3-py 2.3.1 values, 1e-9 |
| M15 | ✅ | `DesignPrimers_RandomTemplate_MatchesPrimer3DesignPrimers`, `DesignPrimers_IndividuallyBestPrimersTmIncompatible_SearchesPairs`, `DesignPrimers_NoPairWithinTmLimit_ReturnsInvalidWithTmMessage`, `DesignPrimers_MatchesBruteForcePrimer3PairSearch` (Differential), `DesignPrimers_KnownTemplate_ProductSizeEqualsSpan` (Properties) | primer3-py design_primers pairs + penalties |
| M16 | ✅ | `PrimerDesigner_PairSearch_Tests.DesignPrimerPairs_Primer3Defaults_MatchesPrimer3Ranks0To4`, `DesignPrimers_ReturnsRankZeroOfDesignPrimerPairs`, `CalculateProductMeltingTemperaturePrimer3_MatchesPrimer3ProductTm` | primer3-py design_primers ranks 0–4 |
| M17 | ✅ | `DesignPrimerPairs_SizeRangesInOrderAndProductTmLimits_MatchesPrimer3`, `DesignPrimerPairs_IncludedRegion_MatchesPrimer3`, `DesignPrimers_DefaultProductSizeRange_Is100To300`, `DesignPrimers_ValidTemplate_ProductWithinDefaultSizeRange`, `DesignPrimers_ProductSizeRange_IsConfigurable`, `DesignPrimers_NoPairWithinTmLimit_ReturnsInvalidWithTmMessage` | primer3-py design_primers |
| M18 | ✅ | `DesignPrimerPairs_NonDefaultPairWeights_MatchesPrimer3PairPenalty` | primer3-py PRIMER_PAIR_k_PENALTY |
| M19 | ✅ | `DesignPrimerPairs_PickInternalOligo_MatchesPrimer3Triples`; MCP `DesignPrimers_PickInternalOligo_MatchesPrimer3` | primer3-py PRIMER_INTERNAL_k_* |
| M20 | ✅ | `DesignPrimers_InvalidOptions_ThrowAsPrimer3DataControl` | primer3-py error strings |
| M21 | ✅ | `PrimerDesigner_EndStabilityWeight_Tests` (PRIMER_WT_END_STABILITY in `EvaluatePrimer` / `DesignPrimerPairs`, thermodynamic + alignment mode, internal oligo; audit round 3, A3-2) | primer3-py PRIMER_LEFT/RIGHT/PAIR_k_PENALTY |
| M23 | ✅ | `PrimerDesigner_ThreePrimeEnd_Tests` (`EvaluatePrimer_GcClamp_MatchesPrimer3` (5), `_MaxEndGc_MatchesPrimer3` (7), `_MaxEndStability_MatchesPrimer3` (6), `Primer3EndChecks_Defaults_AreInactive_AndMatchPrimer3`, `EvaluatePrimer_Avoid3PrimeGC_LibraryRuleKept`, `EvaluatePrimer_IllegalEndChecks_Throw`, `DesignPrimerPairs_GcClamp2_/_MaxEndGc2_/_MaxEndStability_MatchesPrimer3`; audit round 3, A3-7); MCP `EvaluatePrimer_ThreePrimeEndChecks_MatchPrimer3`, `DesignPrimers_ThreePrimeEndChecksAndDistance_MatchPrimer3` | primer3-py 2.3.1 check_primers / design_primers |
| M24 | ✅ | `PrimerDesigner_ThreePrimeEnd_Tests` (`DesignPrimerPairs_MinThreePrimeDistance3_`, `_MinLeft5Right0_`, `_MinThreePrimeDistance0_NoPrimerReused_`, `_DefaultDistance_ReusesPrimers_MatchesPrimer3`, `_MinThreePrimeDistanceBelowMinusOne_Throws`; audit round 3, A3-6); MCP `DesignPrimers_ThreePrimeEndChecksAndDistance_MatchPrimer3` | primer3-py 2.3.1 design_primers ranks 0–4 |
| M22 | ✅ | `PrimerDesigner_ReactionConditions_Tests` (EvaluatePrimer Tm / penalty / ntthal values, thermodynamic + alignment-mode 5-rank designs with product Tm, internal-oligo conditions + PRIMER_INTERNAL_OPT_GC_PERCENT, `DesignProbesPrimer3` internal GC weights, null = defaults, GC optimum inert without weights, 7 illegal-condition cases; audit round 3, A3-1); MCP `DesignPrimers_ReactionConditions_MatchPrimer3`, `EvaluatePrimer_ReactionConditions_MatchPrimer3` | primer3-py 2.3.1 design_primers |
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

| Primer3 alignment mode (audit round 2, A1) | `DesignPrimerPairs_Primer3AlignmentDefaults_MatchesPrimer3`, `_Primer3AlignmentNonDefaultLimitsWeightsAndInternalOligo_MatchesPrimer3`, `_ThermodynamicPerPrimerWeights_MatchesPrimer3`, `DesignPrimers_Primer3AlignmentScreen_ReturnsRankZero`, `DesignPrimerPairs_IllegalComplementarityLimits_Throw`, `DesignProbesPrimer3_AlignmentMode_MatchesPrimer3PickHybProbeOnly` (primer3-py 2.3.1, PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0) | `PrimerDesigner_AlignmentMode_Tests.cs` |

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
- **Structure screens**: heuristic hairpin / primer-dimer / dinucleotide checks (PRIMER-STRUCT-001) instead of Primer3's ntthal Tm limits; no product-size range. In ~30% of random templates Primer3's best pair fails the heuristic hairpin screen, so the chosen pair differs from Primer3's there.
- **Target coordinates**: `targetEnd` is exclusive (target = `[targetStart, targetEnd)`, Primer3 SEQUENCE_TARGET start,length).
