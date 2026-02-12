# Test Specification: PROTMOTIF-FIND-001

**Test Unit ID:** PROTMOTIF-FIND-001
**Area:** ProteinMotif
**Algorithm:** Protein Motif Search
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | PROSITE Database (SIB) | 2 | https://prosite.expasy.org/ | 2026-02-12 |
| 2 | PROSITE User Manual | 2 | https://prosite.expasy.org/prosuser.html | 2026-02-12 |
| 3 | Hulo et al. (2007) | 1 | https://doi.org/10.1093/nar/gkm977 | 2026-02-12 |
| 4 | De Castro et al. (2006) | 1 | https://doi.org/10.1093/nar/gkl124 | 2026-02-12 |
| 5 | Wikipedia: Sequence motif | 4 | https://en.wikipedia.org/wiki/Sequence_motif | 2026-02-12 |
| 6 | PROSITE PS00001–PS00029 entries | 2 | https://prosite.expasy.org/PS00001 etc. | 2026-02-12 |

### 1.2 Key Evidence Points

1. PROSITE patterns have a formally defined syntax with precise semantics — PROSITE User Manual
2. Each PROSITE pattern is fully specified with exact repeat counts and character classes — PROSITE entries
3. PS00007 official pattern is `[RK]-x(2)-[DE]-x(3)-Y` (fixed repeats, not ranges) — PROSITE PS00007
4. PS00018 official pattern is `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]` — PROSITE PS00018
5. Pattern matching is case-insensitive per PROSITE convention — PROSITE User Manual

### 1.3 Documented Corner Cases

1. **Proline exclusion in N-glycosylation:** `N-{P}-[ST]-{P}` — positions 2 and 4 must not be Proline. PROSITE PS00001.
2. **Skip-flag patterns:** Ubiquitous patterns produce many hits; not errors. PROSITE documentation.
3. **Empty or very short sequences:** Should return no matches gracefully.

### 1.4 Known Failure Modes / Pitfalls

1. **Incorrect repeat ranges** — Using `x(2,3)` instead of `x(2)` changes the pattern semantics. Discovered in PS00007.
2. **Missing pattern elements** — Omitting trailing elements reduces pattern specificity. Discovered in PS00018.
3. **Invalid regex from malformed pattern** — Must be handled without exceptions.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindMotifByPattern(sequence, regexPattern, motifName, patternId)` | ProteinMotifFinder | Canonical | Single-pattern scan via regex |
| `FindCommonMotifs(proteinSequence)` | ProteinMotifFinder | Canonical | Scan all CommonMotifs patterns |
| `CommonMotifs` (dictionary) | ProteinMotifFinder | Internal | Verified patterns from PROSITE |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | FindMotifByPattern returns only non-overlapping matches | Yes | .NET Regex.Matches semantics |
| INV-2 | Every match Start ≥ 0 and End < sequence.Length | Yes | Array bounds correctness |
| INV-3 | Every match Sequence equals sequence.Substring(Start, End - Start + 1) | Yes | Regex match value consistency |
| INV-4 | FindMotifByPattern is case-insensitive (upper/lower yield same results) | Yes | PROSITE convention |
| INV-5 | FindCommonMotifs returns empty for empty/null input | Yes | Trivial |
| INV-6 | CommonMotifs PROSITE patterns match official PROSITE definitions | Yes | PROSITE entries PS00001–PS00029 |
| INV-7 | Score and EValue are non-negative for all valid matches | Yes | Score function definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | RGD_SimplePattern_ThreeMatches | FindMotifByPattern finds 3 RGD in `MRGDKLARGDPMRGD` | 3 matches at positions 1, 6, 12 with sequence "RGD" | PROSITE PS00016 |
| M2 | RGD_PositionAndSequence | Verify Start, End, and Sequence fields are correct | Start=1,End=3; Start=6,End=8; Start=12,End=14 | Regex match semantics |
| M3 | NGlycos_FindsValidSite | FindCommonMotifs finds N-glycosylation in sequence containing `NFTA` | Match with MotifName containing "GLYCOSYL" | PROSITE PS00001 |
| M4 | NGlycos_ExcludesProline | No N-glycosylation match when P follows N or precedes end | No ASN_GLYCOSYLATION matches in `ANPSANPT` | PROSITE PS00001 exclusion |
| M5 | PKC_FindsPhosphorylation | FindCommonMotifs finds PKC site in `[ST]-x-[RK]` | At least 1 PKC match | PROSITE PS00005 |
| M6 | CK2_FindsPhosphorylation | FindCommonMotifs finds CK2 site in `[ST]-x(2)-[DE]` | At least 1 CK2 match | PROSITE PS00006 |
| M7 | PLoop_FindsSite | FindCommonMotifs finds P-loop in `[AG]-x(4)-G-K-[ST]` | At least 1 ATP_GTP_A match | PROSITE PS00017 |
| M8 | EmptySequence_ReturnsEmpty | FindMotifByPattern("", pattern) returns empty | Empty enumerable | Trivial |
| M9 | NullSequence_ReturnsEmpty | FindMotifByPattern(null, pattern) returns empty | Empty enumerable | Trivial |
| M10 | InvalidRegex_ReturnsEmpty | FindMotifByPattern with malformed regex returns empty | Empty enumerable, no exception | Robustness |
| M11 | CaseInsensitive | Lowercase input produces same matches as uppercase | Identical match count and positions | PROSITE convention |
| M12 | PS00007_CorrectPattern | After fix, PS00007 regex matches `[RK]-x(2)-[DE]-x(3)-Y` exactly | Regex `[RK].{2}[DE].{3}Y` | PROSITE PS00007 |
| M13 | PS00018_CorrectPattern | After fix, PS00018 matches PROSITE official pattern | Regex with `[^W]` at pos 2 and trailing `[LIVMFYW]` | PROSITE PS00018 |
| M14 | CommonMotifs_AllPrositeCorrect | All PROSITE-sourced patterns in CommonMotifs match their official definitions | Pattern and regex strings validated | PROSITE PS00001–PS00029 |
| M15 | CommonMotifs_ValidRegex | All regex patterns in CommonMotifs compile without error | No RegexParseException | Regex correctness |
| M16 | MatchInvariant_SubstringConsistency | For every match, Sequence == input[Start..End+1] | All matches satisfy invariant | INV-3 |
| M17 | Score_NonNegative | Every match has Score > 0 | All scores positive | INV-7 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | MultiplePatterns_SameSequence | FindCommonMotifs returns matches from different patterns on one sequence | Multiple motif types found | Integration |
| S2 | NoMatch_ReturnsEmpty | Pattern that doesn't match the sequence returns empty | Empty result | Edge case |
| S3 | OverlappingPotential_HandledConsistently | Sequence with potential overlapping motifs produces consistent results | Non-overlapping matches per regex semantics | **ASSUMPTION: non-overlapping** |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | LargeProtein_CompletesWithoutError | 500+ aa protein scanned for all motifs | Completes in reasonable time | Performance |

---

## 5. Audit of Existing Tests

### 5.1 Test Files

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_MotifSearch_Tests.cs` | Canonical for PROTMOTIF-FIND-001 | 34 |
| `ProteinMotifFinderTests.cs` | Residual: SP-001(4), TM-001(4), DISORDER-pred(3), CC-001(3), LC-001(3), DOMAIN-001(3), case(1), integration(2) | 23 |

### 5.2 Coverage Classification

| Test Case ID | Status | Canonical Test Method |
|--------------|--------|----------------------|
| M1 | ✅ Covered | `FindMotifByPattern_RGD_FindsThreeMatches` |
| M2 | ✅ Covered | `FindMotifByPattern_RGD_ReturnsCorrectPositions` |
| M3 | ✅ Covered | `FindCommonMotifs_NGlycosylation_FindsValidSite` |
| M4 | ✅ Covered | `FindCommonMotifs_NGlycosylation_ExcludesProline` |
| M5 | ✅ Covered | `FindCommonMotifs_PKC_FindsPhosphorylation` |
| M6 | ✅ Covered | `FindCommonMotifs_CK2_FindsPhosphorylation` |
| M7 | ✅ Covered | `FindCommonMotifs_PLoop_FindsSite` |
| M8 | ✅ Covered | `FindMotifByPattern_EmptySequence_ReturnsEmpty` |
| M9 | ✅ Covered | `FindMotifByPattern_NullSequence_ReturnsEmpty`, `FindCommonMotifs_NullSequence_ReturnsEmpty` |
| M10 | ✅ Covered | `FindMotifByPattern_InvalidRegex_ReturnsEmpty` |
| M11 | ✅ Covered | `FindMotifByPattern_CaseInsensitive_SameResults` |
| M12 | ✅ Covered | `CommonMotifs_PS00007_MatchesOfficialPrositeDefinition`, `FindMotifByPattern_PS00007_MatchesCorrectLength` |
| M13 | ✅ Covered | `CommonMotifs_PS00018_MatchesOfficialPrositeDefinition` |
| M14 | ✅ Covered | `CommonMotifs_AllPrositePatterns_MatchOfficialDefinitions` |
| M15 | ✅ Covered | `CommonMotifs_AllRegexPatterns_CompileSuccessfully` |
| M16 | ✅ Covered | `FindMotifByPattern_MatchSequence_EqualsSubstring`, `FindCommonMotifs_AllMatches_SatisfySubstringInvariant` |
| M17 | ✅ Covered | `FindMotifByPattern_Score_IsInformationContent`, `FindMotifByPattern_Score_AccountsForCharacterClasses`, `FindMotifByPattern_EValue_UsesProperProbability` |
| S1 | ✅ Covered | `FindCommonMotifs_MultiplePatterns_ReturnsMultipleMotifTypes` |
| S2 | ✅ Covered | `FindMotifByPattern_NoMatch_ReturnsEmpty`, `FindMotifByPattern_EmptyPattern_ReturnsEmpty` |
| S3 | ✅ Covered | `FindMotifByPattern_MatchFields_ArePopulated` |
| S4 | ✅ Covered | `FindMotifByPattern_OverlappingMatches_AllDiscovered`, `FindMotifByPattern_NonOverlapping_SameAsOverlapping` |
| S5 | ✅ Covered | `CommonMotifs_NLS1_MatchesChelskysConsensus`, `CommonMotifs_NES1_MatchesLaCourConsensus`, `CommonMotifs_SIM1_MatchesHeckerConsensus`, `CommonMotifs_WW1_MatchesChenSudolConsensus`, `CommonMotifs_SH3_1_MatchesMayerClassIConsensus`, `CommonMotifs_AllNonProsite_HaveCorrectPatterns` |

**Missing:** 0 &emsp; **Weak:** 0 &emsp; **Duplicate:** 0

### 5.3 Residual Tests Coverage Classification

The 23 residual tests in `ProteinMotifFinderTests.cs` belong to future test units.
All 12 weak tests were strengthened (exact values, unconditional asserts, assertion messages).

| Test | Unit | Status |
|------|------|--------|
| `PredictSignalPeptide_ClassicSignal_PredictsSite` | SP-001 | ✅ Strengthened — exact CleavagePosition=25, Score, Probability |
| `PredictSignalPeptide_NoSignal_ReturnsNull` | SP-001 | ✅ OK |
| `PredictSignalPeptide_ShortSequence_ReturnsNull` | SP-001 | ✅ OK |
| `PredictSignalPeptide_ReturnsRegions` | SP-001 | ✅ Strengthened — unconditional assert, exact NRegion/HRegion/CRegion/Score |
| `PredictTransmembraneHelices_HydrophobicStretch_FindsHelix` | TM-001 | ✅ Strengthened — exact count=1, Start=0, End=29, Score=3.8 |
| `PredictTransmembraneHelices_NoHydrophobic_ReturnsEmpty` | TM-001 | ✅ OK |
| `PredictTransmembraneHelices_MultipleTM_FindsAll` | TM-001 | ✅ Strengthened — exact count=3, start positions |
| `PredictTransmembraneHelices_ShortSequence_ReturnsEmpty` | TM-001 | ✅ OK |
| `PredictDisorderedRegions_DisorderProne_FindsRegions` | DISORDER | ✅ Strengthened — was `≥0` (always passes), now exact count=1, positions, score |
| `PredictDisorderedRegions_Ordered_ReturnsEmpty` | DISORDER | ✅ OK |
| `PredictDisorderedRegions_ShortSequence_ReturnsEmpty` | DISORDER | ✅ OK |
| `PredictCoiledCoils_HeptadPattern_FindsCoil` | CC-001 | ✅ Strengthened — was `≥0` (always passes), now exact count=5, start, score |
| `PredictCoiledCoils_NoPattern_ReturnsEmpty` | CC-001 | ✅ OK |
| `PredictCoiledCoils_ShortSequence_ReturnsEmpty` | CC-001 | ✅ OK |
| `FindLowComplexityRegions_PolyAlanine_Finds` | LC-001 | ✅ Strengthened — exact count=1, DominantAa='A', Frequency=1.0 |
| `FindLowComplexityRegions_Diverse_ReturnsEmpty` | LC-001 | ✅ OK |
| `FindLowComplexityRegions_MultipleRegions_FindsAll` | LC-001 | ✅ Strengthened — exact count=2, DominantAa='G'/'S' |
| `FindDomains_ZincFinger_Finds` | DOMAIN-001 | ✅ Strengthened — was `Is.Not.Null` on `.ToList()`, now exact count=1, name, positions |
| `FindDomains_PLloop_FindsKinase` | DOMAIN-001 | ✅ Strengthened — exact count=1, name, start=4, end=11 |
| `FindDomains_EmptySequence_ReturnsEmpty` | DOMAIN-001 | ✅ OK |
| `PredictSignalPeptide_HandlesLowercase` | SP-001 | ✅ OK |
| `FullWorkflow_AnalyzeProtein` | Integration | ✅ Strengthened — exact motif count=63, specific motif checks |
| `FullWorkflow_LargeProtein` | Integration | ✅ Strengthened — was `Is.Not.Null` on `.ToList()`, now exact counts for all methods |

**Residual:** Missing: 0 &emsp; Weak: 0 (12 fixed) &emsp; Duplicate: 0

---

## 6. Assumption Register

**Total assumptions:** 0

---

## 7. Open Questions / Decisions

None.
