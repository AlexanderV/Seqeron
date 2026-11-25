# Test Specification: SPLICE-DONOR-001

**Test Unit ID:** SPLICE-DONOR-001
**Area:** Splicing
**Algorithm:** Donor (5') Splice Site Detection
**Status:** ✅ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25 (B11 review: Shapiro & Senapathy 1987 score, U12 /RTATCCTY, MaxEntScan mode)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| S1 | Shapiro & Senapathy (1987). Nucleic Acids Res 15(17):7155-7174 | 1 | doi:10.1093/nar/15.17.7155 | 2026-02-12 |
| S2 | Burge, Tuschl & Sharp (1999). The RNA World, CSHL Press | 1 | ISBN 978-0-87969-380-0 | 2026-02-12 |
| S3 | Yeo & Burge (2004). J Comput Biol 11(2-3):377-394 | 1 | doi:10.1089/106652704773135290 | 2026-02-12 |
| S4 | Wikipedia: RNA splicing | 4 | https://en.wikipedia.org/wiki/RNA_splicing | 2026-02-12 |
| S5 | Wikipedia: Spliceosome | 4 | https://en.wikipedia.org/wiki/Spliceosome | 2026-02-12 |

### 1.2 Key Evidence Points

1. Donor site has almost invariant GU dinucleotide at intron start — S1, S4
2. Extended consensus is MAG|GURAGU (positions -3 to +6) — S1
3. Position 0 (G) and +1 (U) are ~100% conserved — S1
4. Position -1 (G) 78%, -2 (A) 58%, -3 C 37% / A 32% (primate, S1 Table 1)
5. GC-AG introns are valid (~0.5-1% of U2-type) — S2
6. U12-type AT-AC introns use AT at donor (~0.3%) — S2
7. S&S score = (t − mint)/(maxt − mint) over the 8-nt subsequence −2..+6 with Table 1 percentages (S1 Methods; worked example AAGTGAGT → 84.3; Table 5 published scores) — reported /100
9. U12-type 5' ss consensus /RTATCCTY (+1..+8), AT-AC and GT-AG U12 introns — Sharp & Burge 1997, S2, Dietrich et al. 1997
8. Higher PWM score indicates stronger splice site — S3

### 1.3 Documented Corner Cases

1. Sequences without any GT/GU dinucleotide — no valid donor sites possible.
2. Sequence shorter than PWM window (~6 nt minimum) — insufficient context for scoring.
3. GC donors are valid but weaker than GT — Burge et al. (1999).
4. Multiple GT occurrences in a sequence — each must be independently evaluated.

### 1.4 Known Failure Modes / Pitfalls

1. Cryptic splice sites: GT dinucleotides in coding exons can score moderately — S1
2. PWM scoring depends critically on context; isolated GT without conserved context should score low — S1, S3

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindDonorSites(sequence, minScore, includeNonCanonical)` | SpliceSitePredictor | **Canonical** | Deep evidence-based testing |
| `FindDonorSites(sequence, minScore, includeNonCanonical, model, category)` | SpliceSitePredictor | **Canonical** | S&S subtable / MaxEntScan model |
| `ScoreDonorSiteShapiroSenapathy`, `ScoreDonorSiteMaxEnt`, `ScoreU12DonorSite` | SpliceSitePredictor | **Internal** | Private; tested via FindDonorSites |
| `ScoreDonorMaxEnt(window)` | SpliceSitePredictor | **Canonical** | Opt-in MaxEntScan score5ss; deep evidence-based testing (Yeo & Burge 2004) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned site has Type == Donor (or U12Donor if non-canonical) | Yes | Algorithm design |
| INV-2 | Score is in [0, 1] range (S&S and U12; MaxEnt = bits) | Yes | (t − mint)/(maxt − mint); matches/8 |
| INV-3 | Confidence is in [0, 1] range | Yes | CalculateConfidence clamps |
| INV-4 | Empty/null input yields no results | Yes | Guard clause |
| INV-5 | A candidate whose scored window runs off the sequence / contains non-ACGTU is not reported (< 8 nt → none) | Yes | S1 scores whole 8-nt subsequences |
| INV-6 | Perfect consensus (CAGGUAAGU) scores ≥ any weaker context | Yes | PWM properties — S1 |
| INV-7 | GC donor score < equivalent GT donor score | Yes | +2 C = 0 % in Table 1 (448/548 vs 1.0) — S1, S2 |
| INV-8 | Motif string is non-empty for every returned site | Yes | GetMotifContext logic |
| INV-9 | All sites are at valid positions within sequence bounds | Yes | Loop bounds |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | CanonicalGT_ConsensusMotif | CAGGUAAGU: every scored position has its top primate percentage | 1 site, pos 3, Score 1.0, Motif CAGGUAAGU | S1 Table 1 |
| M1b | ShapiroSenapathy_ReproducesPublishedTable5Scores | 23 Table 5 subsequences, 6 subtables; worked example AAGTGAGT = 462/548 | score·100 = published (1 dp) | S1 Table 5 / Methods |
| M2 | NoGU_ReturnsEmpty | AAAAACCCCC has no GU → empty result | Empty collection | Trivially correct |
| M3 | EmptyInput_ReturnsEmpty | Empty string → empty result | Empty collection | Trivially correct |
| M4 | SequenceShorterThanScoredWindow | GUAA, AGGUAAG (7 nt) → empty; AGGUAAGU → pos 2 | Empty / 1 site | S1: 8-nt subsequence |
| M5 | StrongVsWeak_ScoreOrdering | CAGGUAAGU 1.0 > UUUGUAAUU 356/548 | exact | S1 |
| M6 | GC_NonCanonical_Detected | CAGGCAAGU with includeNonCanonical=true | 1 Donor, score 448/548 | S1, S2 |
| M7 | GC_NotDetected_WhenCanonicalOnly | CAGGCAAGU with includeNonCanonical=false → no GT present → empty | Empty or no canonical donor | S2 |
| M8 | DNA_T_Equivalence | CAGGTAAGT (T) produces same site count as CAGGUAAGU (U) | Same number of sites | Implementation T→U conversion |
| M9 | LowercaseHandling | cagguaagu lowercase → finds donor site | ≥1 site | Implementation ToUpperInvariant |
| M10 | MultipleSites_AllDetected | 3 GUs | positions {3,7,26}, scores {1, 284/548, 1} | S1 |
| U1 | U12_AT_Donor / U12_GTAG_Donor / U12_PartialMatch | /AUAUCCUU, /GUAUCCUU (also U2 Donor 376/548), /AUAACCAA | U12Donor 1.0, 1.0, 5/8 | Sharp & Burge 1997; Dietrich 1997 |
| X1 | MaxEntScanMode_* | CAGGUAAGU → 10.858313101356437 bits, conf 0.9994616209198077; GC 3.104161570996358; window −3..+6 must fit | exact | S3 (maxentpy score5) |
| ME1 | ScoreDonorMaxEnt_Canonical_10.86 | `ScoreDonorMaxEnt("cagGTAAGT")` == 10.86 bits (2 dp) | 10.86 | S3: maxentpy `score5` docstring (canonical reference) |
| ME2 | ScoreDonorMaxEnt_Canonical_FullPrecision | Full-precision value behind 10.86 | 10.858313 ± 1e-6 | S3: reproduced this session |
| ME3 | ScoreDonorMaxEnt_Stronger_11.08 | `ScoreDonorMaxEnt("gagGTAAGT")` == 11.08 bits | 11.08 / 11.078494 | S3: maxentpy `score5` docstring |
| ME4 | ScoreDonorMaxEnt_Weak_Minus0.12 | `ScoreDonorMaxEnt("taaATAAGT")` == -0.12 bits | -0.12 / -0.116791 | S3: maxentpy `score5` docstring (non-GT) |
| ME5 | ScoreDonorMaxEnt_StrongRanksAboveWeak | strong (10.86) > weak (-0.12) | strong > weak | S3: ordering of documented examples |
| ME6 | ScoreDonorMaxEnt_DnaRnaEquivalence | T-form and U-form windows score identically | equal ± 1e-12 | T==U in rest key / GT model |
| ME7 | ScoreDonorMaxEnt_CaseInsensitive | upper-case == lower-case score | equal ± 1e-12 | Implementation ToUpperInvariant |
| ME8 | ScoreDonorMaxEnt_Null_Throws | null window → ArgumentNullException | throws | Implementation guard |
| ME9 | ScoreDonorMaxEnt_InvalidWindow_Throws | wrong length / non-A/C/G/T(/U) → ArgumentException | throws | Implementation guard (9-nt, alphabet) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | ScoreInRange | All returned sites have Score in [0, 1] | True for every site | INV-2 |
| S2 | ConfidenceInRange | All returned sites have Confidence in [0, 1] | True for every site | INV-3 |
| S3 | MotifNonEmpty | All returned sites have non-empty Motif | True | INV-8 |
| S4 | HighThreshold_ReducesSites | minScore=0.8 returns fewer or equal sites to minScore=0.2 | highCount ≤ lowCount | Score filtering |
| S5 | NullInput_ReturnsEmpty | null input → empty, no exception | Empty collection | Defensive |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | U12_AT_Donor | /AUAUCCUU with includeNonCanonical → U12Donor 1.0 | U12Donor | Sharp & Burge 1997 |
| C2 | GC_Donor_LowerThanGT | GC 448/548 < GT 1.0 — INV-7 | GC score < GT score | S1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs` (18 tests)
- Donor tests fully consolidated from shared `SpliceSitePredictorTests.cs`

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 — Canonical GT detection | ✅ Covered | Exact score=1.0 (9/9), confidence=1.0, position=3, type=Donor |
| M2 — No GU returns empty | ✅ Covered | `FindDonorSites_NoGU_ReturnsEmpty` |
| M3 — Empty input | ✅ Covered | `FindDonorSites_EmptyString_ReturnsEmpty` |
| M4 — Short sequence | ✅ Covered | `FindDonorSites_SequenceShorterThan6_ReturnsEmpty` |
| M5 — Score ordering | ✅ Covered | Exact values: strong=9/9=1.0, weak=5/9; ordering verified |
| M6 — GC non-canonical | ✅ Covered | Exact score=8/9, type=Donor, position=3 |
| M7 — GC not detected canonical-only | ✅ Covered | Clean GC-only sequence (CAGGCAACC), no GU anywhere |
| M8 — DNA T equivalence | ✅ Covered | Count, score=1.0, position equality verified |
| M9 — Lowercase handling | ✅ Covered | Exact count=1, score=1.0, type=Donor |
| M10 — Multiple sites | ✅ Covered | Exact 3 sites at positions {3,7,26} with scores {1.0, 4/9, 1.0} |
| ME1 — MaxEnt canonical 10.86 | ✅ Covered | `ScoreDonorMaxEnt_CanonicalReferenceWindow_Returns10Point86` (round 2 dp == 10.86) |
| ME2 — MaxEnt full precision | ✅ Covered | `ScoreDonorMaxEnt_CanonicalReferenceWindow_MatchesFullPrecision` (10.858313 ± 1e-6) |
| ME3 — MaxEnt 11.08 | ✅ Covered | `ScoreDonorMaxEnt_StrongerSiteReferenceWindow_Returns11Point08` |
| ME4 — MaxEnt -0.12 | ✅ Covered | `ScoreDonorMaxEnt_WeakSiteReferenceWindow_ReturnsMinus0Point12` |
| ME5 — MaxEnt ordering | ✅ Covered | `ScoreDonorMaxEnt_StrongSite_RanksAboveWeakSite` |
| ME6 — DNA/RNA equivalence | ✅ Covered | `ScoreDonorMaxEnt_DnaAndRnaWindows_ProduceIdenticalScores` |
| ME7 — Case insensitive | ✅ Covered | `ScoreDonorMaxEnt_UpperCaseWindow_ScoresIdenticallyToLowerCase` |
| ME8 — Null throws | ✅ Covered | `ScoreDonorMaxEnt_NullWindow_Throws` |
| ME9 — Invalid window throws | ✅ Covered | `ScoreDonorMaxEnt_InvalidWindow_Throws` (8-nt, 10-nt, non-ACGT) |
| S1 — Score range | ✅ Covered | All scores in [0, 1] for multi-site sequence |
| S2 — Confidence range | ✅ Covered | All confidences in [0, 1] |
| S3 — Motif non-empty | ✅ Covered | All motifs non-null and non-empty |
| S4 — Higher threshold | ✅ Covered | minScore=0.8 ≤ minScore=0.2 count |
| S5 — Null input | ✅ Covered | `FindDonorSites_NullInput_ReturnsEmpty` |
| C1 — U12 AT donor | ✅ Covered | AU donor detected as U12Donor type |
| C2 — GC < GT scoring | ✅ Covered | GT=1.0 (9/9), GC=8/9; unconditional assertions |

**Summary:** 0 missing, 0 weak, 0 duplicate. All assertions use theory-derived exact values from MAG|GURAGU consensus.

### 5.3 Consolidation Plan

Consolidation complete. Donor tests removed from `SpliceSitePredictorTests.cs` (retains note at line 11).

### 5.4 Final State

| File | Role | Test Count |
|------|------|------------|
| `SpliceSitePredictor_DonorSite_Tests.cs` | Canonical for SPLICE-DONOR-001 | 18 |
| `SpliceSitePredictorTests.cs` | Remaining tests for SPLICE-ACCEPTOR/PREDICT | Reduced |

---

## 6. Assumption Register

**Total assumptions:** 0

- The binary IUPAC consensus fraction used before 2026-11-25 (unsourced as a score) was replaced by the published
  Shapiro & Senapathy (1987) score (B11 F23). The U12 score is the /RTATCCTY consensus match fraction (no
  bundleable trained U12 PWM — see Evidence). Default `minScore` 0.5 = S&S score 50 (caller-chosen threshold).

## 7. Open Questions / Decisions

None. Evidence is sufficient to define correct behavior for all MUST tests.
