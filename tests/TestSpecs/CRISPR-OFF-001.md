# TestSpec: CRISPR-OFF-001 - Off-Target Analysis

## Test Unit Identification
- **Test Unit ID**: CRISPR-OFF-001
- **Algorithm Group**: MolTools
- **Algorithm Name**: Off_Target_Analysis
- **Canonical Methods**:
  - `CrisprDesigner.FindOffTargets(string, DnaSequence, int, CrisprSystemType)`
  - `CrisprDesigner.CalculateSpecificityScore(string, DnaSequence, CrisprSystemType)`
- **Also owned by this unit** (published off-target models in the same class):
  - `CrisprDesigner.CalculateMitHitScore(string, string)` — MIT/Hsu 2013 single hit
  - `CrisprDesigner.CalculateSaCas9HitScore(string, string)` — Tycko et al. 2018 single hit
  - `CrisprDesigner.CalculateCfdScore(string, string, string)` — Doench 2016 CFD single hit
  - `CrisprDesigner.CalculateMitSpecificityScore(IEnumerable<double>)` and `(string, DnaSequence, int, CrisprSystemType)`
  - `CrisprDesigner.CalculateCfdSpecificityScore(IEnumerable<double>, CfdAggregation)` and `(string, DnaSequence, int, CrisprSystemType, CfdAggregation)`
  - `CfdAggregation`
- **Related Documentation**: [Off_Target_Analysis.md](../docs/algorithms/MolTools/Off_Target_Analysis.md),
  [CRISPR-OFF-001-Evidence.md](../../docs/Evidence/CRISPR-OFF-001-Evidence.md)
- **Last reviewed**: 2026-10-14 (campaign `docs/Validation/review-2026-09`, batch B08)

---

## Evidence Summary

### Primary Sources
| Source | Key Information | URL/Reference |
|--------|-----------------|---------------|
| Wikipedia: Off-target genome editing | Off-target mutations occur with 3-5 mismatches; seed sequence (10-12nt adjacent to PAM) critical for specificity; PAM-distal mismatches more tolerated | https://en.wikipedia.org/wiki/Off-target_genome_editing |
| Hsu et al. (2013), Nature Biotechnology | PAM-proximal 8-12bp defines specificity; single-base specificity ranges 8-14bp; 2+ mismatches in PAM-proximal region reduce activity; 3+ interspaced or 5+ concatenated mismatches eliminate cleavage | PMC3969858 |
| Fu et al. (2013), Nature Biotechnology | >50% of CRISPR mutations can occur off-target; truncated gRNAs (17-18bp) increase specificity | PMID 23792628 |
| Doench et al. (2016), Nature Biotechnology | CFD (Cutting Frequency Determination): per-position mismatch percent-activity matrix × PAM-activity table | PMID 26780180 |
| Tycko J, Wainberg M, Marinov GK, et al. (2018), Nature Communications 9:2962 | SaCas9-specific trained per-position mismatch coefficients (pairwise library screen) | doi:10.1038/s41467-018-05391-2 |
| CRISPOR `crispor.py` (master) | Reference implementation of every model and both guide-level aggregates (`calcHitScore`, `calcSaHitScore`, `calcCfdScore`, `calcMitGuideScore`, `calcCfdGuideScore`, `annotateOfftargets`) | https://raw.githubusercontent.com/maximilianh/crisporWebsite/master/crispor.py |

> **Provenance note (WP-SRC, 2026-10-14).** The **Wikipedia: Off-target genome editing** row is on a
> host denied by the organization's egress policy (`en.wikipedia.org` — agent proxy answers 403 to
> CONNECT); **no unit of this campaign has opened it**, so it is **unverifiable in this environment**
> and must not be read as checked. It is retained, not deleted: every parameter it is cited for is
> also carried by a reachable source in the same table — the seed-region and mismatch-tolerance
> claims by Hsu et al. 2013, and all scoring models by the executed CRISPOR reference. See
> [`docs/Evidence/B08-Source-Provenance.md`](../../docs/Evidence/B08-Source-Provenance.md) §1, §7c.

### Evidence-Backed Parameters
| Parameter | Evidence Value | Implementation Value | Status |
|-----------|----------------|---------------------|--------|
| Max mismatches for off-target | 3-5 mismatches | 0-5 (configurable) | ✅ Aligned |
| Seed region length (Cas9) | 8-12bp PAM-proximal (Hsu 2013) | 12bp | ✅ Conservative upper bound |
| Off-target activity threshold | Mismatches >0 | Mismatches >0 | ✅ Aligned |
| PAM requirement | NGG or NAG for off-targets | PAM required | ✅ Aligned |

### Unsourced element (declared, not fixed)

`CalculateSpecificityScore` = `max(0, 100 − Σ(5 per seed mismatch + 2 per distal mismatch))`.
Hsu 2013 motivates only the *direction* of the seed weighting; **no source was found for the weights
`5`/`2`, the 12-nt seed width, the linear sum or the subtraction from 100**. Per the campaign's rule
for unsourced convenience metrics it is kept and declared honestly (XML doc, algorithm doc,
`crispr_specificity_score` MCP description, Evidence file), and the published aggregates in the same
class are offered alongside it everywhere, including over MCP as `crispr_mit_specificity_score` and
`crispr_cfd_specificity_score`. Test `CalculateSpecificityScore_HeuristicAndMitAggregate_AreDistinct`
pins that the two produce visibly different numbers on the same input.

### Key Insights from Evidence
1. **Off-targets require mismatches**: By definition, off-target sites have sequence differences from the guide (Hsu 2013)
2. **Seed region importance**: Mismatches in the seed region (PAM-proximal) are less tolerated but still contribute to off-target scoring (Wikipedia, Hsu 2013)
3. **Position-dependent scoring**: Mismatches at different positions have different impacts on cleavage activity (Hsu 2013)
4. **PAM tolerance**: Off-targets can occur at NAG sites as well as NGG (Hsu 2013: "5× lower efficiency at NAG")

---

## Test Categories

### MUST Tests (Critical Functionality)

#### M-001: FindOffTargets - Empty Guide Throws Exception
**Evidence**: Defensive programming; null/empty input undefined
**Input**: Empty string "", valid genome, maxMismatches=3
**Expected**: Throws `ArgumentNullException`
**Existing Test**: `FindOffTargets_EmptyGuide_ThrowsException` ✅

#### M-002: FindOffTargets - Null Genome Throws Exception
**Evidence**: Defensive programming
**Input**: Valid guide, null genome, maxMismatches=3
**Expected**: Throws `ArgumentNullException`
**Existing Test**: `FindOffTargets_NullGenome_ThrowsException` ✅

#### M-003: FindOffTargets - Invalid MaxMismatches Throws Exception
**Evidence**: Per Hsu et al., practical limit is 5 mismatches for detectable off-targets
**Input**: Valid guide, valid genome, maxMismatches > 5 (e.g., 10)
**Expected**: Throws `ArgumentOutOfRangeException`
**Existing Test**: `FindOffTargets_InvalidMaxMismatches_ThrowsArgumentOutOfRangeException` ✅

#### M-003b: FindOffTargets - Guide Length Mismatch Throws Exception
**Evidence**: Hsu 2013 — SpCas9 uses 20bp guide; Cas12a uses 23bp. Mismatched lengths are invalid.
**Input**: 12bp guide with SpCas9, 24bp guide with SpCas9, 20bp guide with Cas12a
**Expected**: Throws `ArgumentException`
**Existing Test**: `FindOffTargets_GuideLengthMismatch_ThrowsArgumentException` ✅

#### M-004: FindOffTargets - Exact Match Returns Empty (Not Off-Target)
**Evidence**: Hsu 2013 - off-targets are sites with mismatches, exact matches are on-targets
**Input**: Guide with exact match in genome with PAM
**Expected**: No off-target returned for exact match (mismatches > 0 required)
**Existing Test**: `FindOffTargets_ExactMatch_NotReturnedAsOffTarget` ✅ (asserts `Is.Empty`)

#### M-005: FindOffTargets - Single Mismatch Within Max Returns Off-Target
**Evidence**: Hsu 2013 - single mismatches tolerated, especially in PAM-distal region
**Input**: Guide with 1 mismatch site in genome
**Expected**: Off-target returned with Mismatches=1
**Existing Test**: ✅ Covered

#### M-006: FindOffTargets - Max Mismatches Respected
**Evidence**: Algorithm specification
**Input**: Guide, genome with sites having various mismatch counts, maxMismatches=2
**Expected**: All returned off-targets have Mismatches ≤ 2
**Existing Test**: `FindOffTargets_MaxMismatchesRespected` ✅

#### M-007: FindOffTargets - MismatchPositions Reported Correctly
**Evidence**: Hsu 2013 - position of mismatches affects activity; must be tracked
**Input**: Guide with known mismatch site
**Expected**: MismatchPositions.Count == Mismatches; positions correct
**Existing Test**: `FindOffTargets_MismatchPositions_CountMatchesMismatches` ✅ (exact positions [0,1] and score=4.0)

#### M-008: FindOffTargets - Requires PAM at Off-Target Site
**Evidence**: Wikipedia, Hsu 2013 - PAM is required for cleavage
**Input**: Genome with similar sequence but no PAM
**Expected**: No off-target returned
**Existing Test**: ✅ Covered

#### M-009: FindOffTargets - Searches Both Strands
**Evidence**: CRISPR can target either strand
**Input**: Guide with off-target on reverse strand only
**Expected**: Off-target returned with IsForwardStrand=false
**Existing Test**: ✅ Covered

#### M-010: CalculateSpecificityScore - Returns Value 0-100
**Evidence**: Score should be normalized percentage
**Input**: Any valid guide and genome
**Expected**: 0 ≤ score ≤ 100
**Invariant Test**: Both existing tests verify bounds ✅

#### M-011: CalculateSpecificityScore - No Off-Targets Returns High Score
**Evidence**: Few off-targets = high specificity
**Input**: Guide with no off-targets in small genome
**Expected**: Score = 100 (or near-maximum)
**Existing Test**: `CalculateSpecificityScore_NoOffTargets_ReturnsHigh` ✅

#### M-012: CalculateSpecificityScore - More Off-Targets Reduces Score
**Evidence**: More off-targets = lower specificity
**Input**: Guide with multiple off-target sites
**Expected**: Score < 100
**Existing Test**: `CalculateSpecificityScore_ManyOffTargets_ReturnsLower` ✅

---

### SHOULD Tests (Important Quality)

#### S-001: FindOffTargets - Seed Region Mismatches Scored Higher
**Evidence**: Hsu 2013 - PAM-proximal (seed) mismatches less tolerated but reduce off-target activity more
**Input**: Guide with off-target having seed region mismatch vs distal mismatch
**Expected**: Seed mismatch has higher OffTargetScore (more penalty)
**Existing Test**: ✅ Covered

#### S-002: FindOffTargets - Multiple Mismatches Summed
**Evidence**: Hsu 2013 - aggregate effect of multiple mismatches
**Input**: Guide with off-target having 2-3 mismatches
**Expected**: Mismatches count and positions all reported
**Existing Test**: ✅ Covered

#### S-003: FindOffTargets - Cas12a System Supported
**Evidence**: Different CRISPR systems have different characteristics
**Input**: Guide, genome, Cas12a system type
**Expected**: Off-targets found using Cas12a PAM (TTTV) and guide length (23bp)
**Existing Test**: ✅ Covered

#### S-004: CalculateSpecificityScore - Seed Mismatches Penalized More
**Evidence**: Implementation uses position-dependent scoring
**Input**: Controlled off-targets with seed vs non-seed mismatches
**Expected**: Seed mismatches contribute more to penalty
**Existing Test**: `CalculateSpecificityScore_SeedMismatch_LowerThanDistal` ✅ (exact scores: distal=98.0, seed=95.0)

---

### COULD Tests (Nice to Have)

#### C-001: FindOffTargets - NAG PAM Sites Detected
**Evidence**: Hsu 2013 - NAG PAM also causes off-target activity (5× lower)
**Input**: Genome with NAG PAM sites
**Expected**: Off-targets detected (implementation-dependent)
**Note**: Current implementation uses standard PAM, NAG support via SpCas9_NAG system

#### C-002: FindOffTargets - Performance Reasonable for Moderate Genomes
**Evidence**: O(n × m) complexity stated
**Input**: ~10kb genome, 20bp guide
**Expected**: Returns in reasonable time (<1s)

---

## Invariants

1. **Off-target definition**: All returned sites have Mismatches > 0 (exact matches are on-targets)
2. **Mismatch bound**: All returned sites have Mismatches ≤ maxMismatches parameter
3. **Mismatch position consistency**: MismatchPositions.Count == Mismatches
4. **PAM requirement**: All off-targets are at valid PAM sites
5. **Score bounds**: 0 ≤ SpecificityScore ≤ 100
6. **Determinism**: Same inputs produce same outputs

---

## Audit of Existing Tests

### Current State (CrisprDesigner_OffTarget_Tests.cs)

| Test | Coverage | Quality | Assertions |
|------|----------|---------|------------|
| `FindOffTargets_EmptyGuide_ThrowsArgumentNullException` | M-001 | ✅ Strong | Throws `ArgumentNullException` |
| `FindOffTargets_NullGenome_ThrowsArgumentNullException` | M-002 | ✅ Strong | Throws `ArgumentNullException` |
| `FindOffTargets_InvalidMaxMismatches_ThrowsArgumentOutOfRangeException` | M-003 | ✅ Strong | 3 TestCases: -1, 6, 10 |
| `FindOffTargets_GuideLengthMismatch_ThrowsArgumentException` | M-003b | ✅ Strong | 3 TestCases: short/long SpCas9, wrong-length Cas12a |
| `FindOffTargets_ExactMatch_NotReturnedAsOffTarget` | M-004 | ✅ Strong | `Is.Empty` (exact match filtered out) |
| `FindOffTargets_SingleMismatch_ReturnsOffTarget` | M-005 | ✅ Strong | Count=1, position=[0], strand=forward, score=2.0 |
| `FindOffTargets_MaxMismatchesRespected_AllResultsWithinLimit` | M-006 | ✅ Strong | All ≤ maxMismatches |
| `FindOffTargets_MismatchPositions_CountMatchesMismatches` | M-007 | ✅ Strong | Count=1, positions=[0,1], score=4.0 |
| `FindOffTargets_MismatchPositions_ContainsCorrectPositions` | M-007b | ✅ Strong | Count=1, positions=[0] |
| `FindOffTargets_NoPam_NoOffTargetReturned` | M-008 | ✅ Strong | `Is.Empty` |
| `FindOffTargets_ReverseStrand_ReturnsOffTarget` | M-009 | ✅ Strong | Count=1, reverse strand, mismatches=1, position=[19], score=5.0 |
| `CalculateSpecificityScore_ReturnsValueInValidRange` | M-010 | ✅ Strong | 0 ≤ score ≤ 100 |
| `CalculateSpecificityScore_NoOffTargets_Returns100` | M-011 | ✅ Strong | Exact value 100 |
| `CalculateSpecificityScore_WithOffTargets_ScoreReducedFromMaximum` | M-012 | ✅ Strong | Exact value 98.0 (1 distal mismatch, penalty=2) |
| `FindOffTargets_SeedMismatch_HigherOffTargetScore` | S-001 | ✅ Strong | Distal=2.0, Seed=5.0, seed > distal |
| `CalculateSpecificityScore_SeedMismatch_LowerThanDistal` | S-004 | ✅ Strong | Distal score=98.0, Seed score=95.0, seed < distal |
| `FindOffTargets_MultipleMismatches_AllReported` | S-002 | ✅ Strong | Count=1, positions=[0,1,2], score=6.0 |
| `FindOffTargets_Cas12a_UsesCorrectParameters` | S-003 | ✅ Strong | Count=1, mismatches=1, positions=[0], score=5.0 |
| `FindOffTargets_EmptyGenome_ReturnsEmpty` | Edge | ✅ Strong | `Is.Empty` |
| `FindOffTargets_GenomeTooShort_ReturnsEmpty` | Edge | ✅ Strong | `Is.Empty` |
| `FindOffTargets_MaxMismatchesZero_ReturnsEmpty` | Edge | ✅ Strong | `Is.Empty` |
| `FindOffTargets_OffTargetScore_IsNonNegative` | Invariant | ✅ Strong | All scores ≥ 0 |
| `FindOffTargets_SameInput_DeterministicOutput` | Invariant | ✅ Strong | Two runs produce identical results |

**Total: 24 tests (26 test executions with parameterized cases)**

---


### Published-model tests added 2026-10-14 (B08, CRISPR-OFF-001)

Every expected value below comes from an independent Python run of a verbatim transcription of the
reference function, never from the C# code. Full tables of opened sources and reference values:
[CRISPR-OFF-001-Evidence.md](../../docs/Evidence/CRISPR-OFF-001-Evidence.md).

#### F9 regression — the MIT aggregate silently returned 100 for non-20-nt systems

| Test | Locks |
|---|---|
| `CalculateMitSpecificityScore_SaCas9_21ntGuide_RealOffTarget_IsNot100` | SaCas9 21-nt guide, one off-target → `95.91285612347565` (not 100) |
| `CalculateMitSpecificityScore_SaCas9_21ntGuide_TwoMismatches_MatchesReference` | two mismatches → `97.54818977171612` |
| `CalculateMitSpecificityScore_SaCas9_MismatchAtFivePrimeEnd_ScoredAsPerfectHit` | 5'-most mismatch contributes nothing → hit 100, aggregate `50.0` |
| `CalculateMitSpecificityScore_Cas12aClass_WithRealOffTargets_Throws` | real off-targets exist **and** both published aggregates throw |
| `CalculateMitSpecificityScore_PamBeforeTargetSystems_AlwaysThrow` | Cas12a, AsCas12a, LbCas12a, CasX all refused |
| `CalculateSpecificityScore_HeuristicAndMitAggregate_AreDistinct` | heuristic `98.0` vs MIT `50.0` on one input |

#### F14 — SaCas9 model (Tycko et al. 2018)

| Test | Locks |
|---|---|
| `CalculateSaCas9HitScore_PerfectMatch_Returns100` | the reference script's first assert |
| `CalculateSaCas9HitScore_ReferenceOracles_Match` | `4.1293256150042446`, `13.4419832027506` (the reference's own asserted values) |
| `CalculateSaCas9HitScore_CoefficientOrientation_NotReversed` | `dT:rC` at index 1 = `54.6734106002513` vs index 20 = `10.569373315278801` |
| `CalculateSaCas9HitScore_MismatchAtFivePrimeMostPosition_Ignored` | index 0 skipped |
| `CalculateSaCas9HitScore_InvalidInputs_Throw` | non-21-nt, gap `-`, non-ACGT, null/empty all refused |
| `CalculateSaCas9HitScore_LowerCaseInput_TreatedAsUpperCase` | case folding |

#### CFD aggregation (both published forms)

| Test | Locks |
|---|---|
| `CfdSpecificity_MitStyle_TwoHits_MatchesReference` | raw `{0.5, 0.25}` → `99.25558312655087` |
| `CfdSpecificity_ParkinsonNormalized_TwoHits_MatchesReference` | raw `{0.5, 0.25}` → `133.33333333333334` (i.e. `100/Σ`, not `100×100/Σ`) |
| `CfdSpecificity_SingleUnitHit_BothAggregationsMatchReference` | raw `{1.0}` → `99.00990099009901` / `100.0` |
| `CfdSpecificity_MitStyle_EmptyAndZeroSum_Returns100` | Σ = 0 is defined for `MitStyle` |
| `CfdSpecificity_ParkinsonNormalized_ZeroSum_ThrowsArgumentOutOfRange` | Σ = 0 is undefined for `ParkinsonNormalized` |
| `CfdSpecificity_Genome_ReverseStrandHit_UsesProtospacerStrandPam` | reverse-strand hit with PAM `AAG` → `99.74141115651946` / `385.71428610000004` |
| `CfdSpecificity_Genome_SaCas9_Refused` | CFD needs a 20-nt protospacer |
| `CfdSpecificity_Genome_MaxMismatchesOutOfRange_ThrowsArgumentOutOfRange` | F10 regression: `-1`, `6`, `99` all throw in both genome aggregates |
| `CfdSpecificity_Genome_NoOffTargets_Returns100`, `CfdSpecificity_NullScores_ThrowsArgumentNullException`, `CfdSpecificity_Genome_InvalidArguments_Throw` | boundaries and validation |

#### `FindOffTargets` consistency and ordering

| Test | Locks |
|---|---|
| `FindOffTargets_MixedCaseInput_CountPositionsAndPenaltyAgree` | INV-04 for a lower-case guide + genome (`Mismatches` = `MismatchPositions.Count` = 1, penalty `2.0`) |
| `FindOffTargets_EnumerationOrder_ForwardSitesThenReverseSites` | INV-05 with one hit on each strand |
| `DnaSequence_RejectsAmbiguityCodes_SoNonAcgtOffTargetsCannotReachTheScan` | why the reference's `[^ATCG]` exclusion is a defensive mirror here |
| `FindOffTargets_Mismatches_EqualCanonicalHammingDistance` | dedup lock: the reported `Mismatches` equals canonical `SequenceExtensions.HammingDistance` on the same pair, for mixed-case input |
| `CalculateCfdScore_NonAcgtBases_RefusedAlthoughCanonicalComplementAcceptsThem` | dedup lock, other direction: `U`/`N`/`R` are refused by CFD although canonical `GetComplementBase` maps them (`U`→`A`, `N`→`N`, `R`→`Y`), so the strict `CfdComplement` cannot be "cleaned up" into a silently wrong score |
| `FindOffTargets_TwoPerfectMatches_BothDroppedGenomeWide` | the recorded on-target divergence: two perfect matches → no off-targets and both published aggregates return 100, where CRISPOR would keep the second and report 50 |

#### Updated invariants

4. **Mismatch-position consistency** — `MismatchPositions.Count == Mismatches`, and the off-target
   penalty is computed from the same mismatch set, so the three can never disagree.
6. **Determinism** — forward-strand hits in ascending position, then reverse-strand hits.
7. **Published aggregates refuse what the reference does not score** — PAM-before-target systems
   (MIT and CFD) and non-20-nt protospacers (CFD) throw rather than returning a value.

#### Recorded divergences pinned by tests (full list in the Evidence file)

- **On-target skip is wider here.** CRISPOR's skip (`crispor.py` L1920-1929) drops only the 0-mismatch
  hits inside a caller-supplied on-target interval — or, with none, the first 0MM hit it meets — and
  scores any further perfect match as a full `calcHitScore = 100` off-target. This API takes a bare
  reference sequence and no locus, so it drops every exact match genome-wide; a guide with a duplicated
  perfect site therefore scores *more* specific here. Pinned by
  `FindOffTargets_TwoPerfectMatches_BothDroppedGenomeWide`.
- **CFD refuses SaCas9 where the reference returns a vacuous 100.** The reference's per-hit
  `cfdScore = -1` for SaCas9 (L1943-1945) means every hit is dropped by its own `if cfdScore != -1`
  filter (L1961-1962), so its `guideCfdScore` is `calcMitGuideScore(sum([])) = 100` — not a sentinel
  (only Cpf1 gets `guideCfdScore = -1`, L1986-1988). Returning 100 for a guide whose off-targets were
  never CFD-scored is the F9 failure mode, so this library throws. Pinned by
  `CfdSpecificity_Genome_SaCas9_Refused`.

### MCP surface

`crispr_specificity_score` (heuristic, description now says so), plus the two tools added in this
unit: `crispr_mit_specificity_score`, `crispr_cfd_specificity_score`. Tests:
`tests/Seqeron/Seqeron.Mcp.MolTools.Tests/CrisprPublishedSpecificityScoreTests.cs`.

---

## Definition of Done

- [x] Evidence documented from authoritative sources
- [x] All MUST tests implemented and passing
- [x] Test file follows naming convention
- [x] No warnings in test project
- [x] Invariants verified through tests
- [x] Processing Registry updated
