# Test Specification: SPLICE-ACCEPTOR-001

**Test Unit ID:** SPLICE-ACCEPTOR-001
**Area:** Splicing
**Algorithm:** Acceptor Site Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-16

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Shapiro & Senapathy (1987). Nucleic Acids Res 15(17):7155–7174 | 1 | — | 2026-02-12 |
| 2 | Burge, Tuschl & Sharp (1999). The RNA World, 2nd ed. CSHL Press | 1 | — | 2026-02-12 |
| 3 | Yeo & Burge (2004). J Comput Biol 11(2–3):377–394 | 1 | — | 2026-02-12 |
| 4 | Patel & Steitz (2003). Nat Rev Mol Cell Biol 4(12):960–970 | 1 | — | 2026-02-12 |
| 5 | Hall & Padgett (1994). J Mol Biol 239(3):357–365 | 1 | — | 2026-03-16 |
| 6 | Jackson (1991). Nucleic Acids Res 19(14):3795–3798 | 1 | — | 2026-03-16 |
| 7 | Dietrich, Incorvaia & Padgett (1997). Molecular Cell 1(1):151–160 | 1 | — | 2026-03-16 |
| 8 | Wikipedia — RNA splicing | 4 | https://en.wikipedia.org/wiki/RNA_splicing | 2026-02-12 |
| 9 | Wikipedia — Polypyrimidine tract | 4 | https://en.wikipedia.org/wiki/Polypyrimidine_tract | 2026-02-12 |
| 10 | Wikipedia — Minor spliceosome | 4 | https://en.wikipedia.org/wiki/Minor_spliceosome | 2026-03-16 |

### 1.2 Key Evidence Points

1. 3' splice site consensus is (Y)nNCAG|G — AG dinucleotide is almost invariant — Shapiro & Senapathy (1987)
2. PPT is 15–20 nt pyrimidine-rich region upstream of AG — Lodish et al. (2004) via Wikipedia
3. PPT quality (continuous pyrimidines) determines splice site strength — Burge et al. (1999)
4. U12-type introns use AC instead of AG at 3' splice site, with YCCAC consensus — Hall & Padgett (1994), Patel & Steitz (2003)
5. Position -2 = A (100%), position -1 = G (100%), position -3 = C (~70%) — Shapiro & Senapathy (1987)
6. First exonic nucleotide (position 0) favors G (~50%) — Shapiro & Senapathy (1987)
7. Terminal dinucleotides alone do not distinguish U2- from U12-type introns — Dietrich et al. (1997)

### 1.3 Documented Corner Cases

1. **U12 YCCAC acceptor**: Minor spliceosome uses AC instead of AG, with YCCAC 3' consensus — Hall & Padgett (1994), Patel & Steitz (2003)
2. **Weak PPT**: Interrupted PPT reduces score; may be skipped in alternative splicing — Burge et al. (1999)
3. **Cryptic AG sites**: AG dinucleotides in intronic sequence can mimic acceptor sites; context and PPT quality distinguish real from cryptic — Shapiro & Senapathy (1987)

### 1.4 Known Failure Modes / Pitfalls

1. **Insufficient upstream context**: Acceptor scoring requires PPT assessment; sequences < 20 nt cannot provide this — implementation guard
2. **T/U ambiguity**: DNA (T) and RNA (U) representations must be handled equivalently — implementation converts T→U

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindAcceptorSites(sequence, minScore, includeNonCanonical)` | SpliceSitePredictor | Canonical | Main scanning method |
| `ScoreAcceptorSite(sequence, position)` | SpliceSitePredictor | Internal | Private; tested via FindAcceptorSites |
| `ScoreU12AcceptorSite(sequence, position)` | SpliceSitePredictor | Internal | Private; tested via FindAcceptorSites with includeNonCanonical=true |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Score ∈ [0, 1] for all returned sites | Yes | Normalization formula clamps to [0, 1] |
| INV-2 | Confidence ∈ [0, 1] for all returned sites | Yes | CalculateConfidence clamps to [0, 1] |
| INV-3 | All canonical sites have Type = Acceptor | Yes | Implementation assigns SpliceSiteType.Acceptor for AG |
| INV-4 | All U12 sites have Type = U12Acceptor | Yes | Implementation assigns SpliceSiteType.U12Acceptor for AC |
| INV-5 | Position = i + 1 (index of G in AG — last intronic nucleotide) | Yes | Implementation: `Position: i + 1` |
| INV-6 | Motif is non-empty for all returned sites | Yes | GetMotifContext extracts surrounding context |
| INV-7 | Empty/null input → empty result | Yes | Guard: `if (string.IsNullOrEmpty(sequence) \|\| sequence.Length < 20) yield break` |
| INV-8 | Higher minScore → subset of lower minScore results | Yes | Score filtering: `if (score >= minScore)` |
| INV-9 | Strong PPT score ≥ weak PPT score for same AG context | Yes | PPT pyrimidine fraction contributes to score — Burge et al. (1999) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Canonical AG detected | Sequence with strong PPT followed by CAG — exactly 1 site, Type=Acceptor, position=18, score≈0.84 | 1 site, Type=Acceptor, pos=18, score≈0.839 | Shapiro & Senapathy (1987) — AG consensus |
| M2 | No AG returns empty | All-pyrimidine sequence (no AG) — no acceptor sites | Empty result | Trivially correct |
| M3 | Empty input returns empty | `""` and `null` inputs should return empty | Empty result | Guard behavior |
| M4 | Short sequence returns empty | Sequence < 20 chars | Empty result | Guard: length < 20 |
| M5 | Strong PPT > weak PPT | Continuous pyrimidine PPT scores higher (>0.7) than purine-interrupted context (<0.7) | Strong score > weak score, both found | Burge et al. (1999) — PPT quality |
| M6 | Score range [0, 1] | All returned scores within [0, 1] | Scores ∈ [0, 1] | INV-1 — normalization |
| M7 | Confidence range [0, 1] | All returned confidences within [0, 1] | Confidence ∈ [0, 1] | INV-2 — CalculateConfidence |
| M8 | DNA T equivalence | DNA input (T instead of U) produces same results | Same site count and scores | Implementation T→U conversion |
| M9 | Case insensitivity | Lowercase input produces same results as uppercase | Same site count, scores, and positions | Implementation ToUpperInvariant |
| M10 | Multiple AG sites | Sequence with two AG dinucleotides in scannable range → exactly 2 detected at correct positions | 2 sites, positions verified, first scores higher | Scanning algorithm |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | U12 YCCAC detected | includeNonCanonical=true with YCCAC consensus → exactly 1 U12Acceptor, position=17, score=1.0 | 1 U12Acceptor, pos=17, score=1.0 | Hall & Padgett (1994), Patel & Steitz (2003) |
| S2 | U12 YCCAC excluded | includeNonCanonical=false (default) with YCCAC → no U12Acceptor sites | No U12Acceptor in result | Default parameter behavior |
| S3 | Position of G in AG | Returned position = i + 1 = index of G (last intronic nucleotide) | Position = 16 for AG at index 15 | INV-5 |
| S4 | Motif non-empty | All returned sites have non-empty Motif containing AG context | Non-empty motif | INV-6 |
| S5 | Threshold filtering | minScore=0.8 returns subset of minScore=0.2 results | High-threshold count ≤ low-threshold count | INV-8 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | AG before position 15 not detected | AG at position < 15 in sequence not found | Not in results | Implementation scan starts at i=15 |
| C2 | PPT pyrimidine count correlates | Independent PPT scoring helper matches expected trend | Score reflects PPT quality | PPT quality principle — Burge et al. (1999) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs` — 17 evidence-based tests
- **Shared file:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictorTests.cs` — integration tests (belong to SPLICE-PREDICT-001)
- **Old acceptor tests:** removed from shared file (replaced by canonical tests)

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (canonical AG) | ✅ | Exact count, position, score, type — `FindAcceptorSites_CanonicalAG_WithStrongPPT_FindsAcceptorSite` |
| M2 (no AG → empty) | ✅ | `FindAcceptorSites_NoAGDinucleotide_ReturnsEmpty` |
| M3 (empty/null) | ✅ | `FindAcceptorSites_EmptyInput_ReturnsEmpty` |
| M4 (short sequence) | ✅ | `FindAcceptorSites_ShortSequence_ReturnsEmpty` |
| M5 (strong vs weak PPT) | ✅ | Both sites found unconditionally, score ranges verified — `FindAcceptorSites_StrongPPT_ScoresHigherThanWeakContext` |
| M6 (score range) | ✅ | `FindAcceptorSites_AllScores_InZeroOneRange` |
| M7 (confidence range) | ✅ | `FindAcceptorSites_AllConfidences_InZeroOneRange` |
| M8 (DNA T equivalence) | ✅ | `FindAcceptorSites_DnaInput_ProducesSameResults` |
| M9 (case insensitivity) | ✅ | Count, score, and position equality — `FindAcceptorSites_LowercaseInput_ProducesSameResults` |
| M10 (multiple AG sites) | ✅ | Exact count=2, positions and score ordering — `FindAcceptorSites_MultipleAGSites_FindsAll` |
| S1 (U12 YCCAC detected) | ✅ | Exact count=1, position, score=1.0 — `FindAcceptorSites_U12AcceptorYCCAC_DetectedWhenNonCanonicalEnabled` |
| S2 (U12 YCCAC excluded) | ✅ | `FindAcceptorSites_U12AcceptorYCCAC_ExcludedByDefault` |
| S3 (G position in AG) | ✅ | Exact count=1, position=16, unconditional — `FindAcceptorSites_ReturnsPositionAfterAG` |
| S4 (motif non-empty) | ✅ | `FindAcceptorSites_MotifIsNonEmpty` |
| S5 (threshold filtering) | ✅ | `FindAcceptorSites_HigherThreshold_ProducesSubset` |
| C1 (AG before scan start) | ✅ | `FindAcceptorSites_AGBeforeScanStart_NotDetected` |
| C2 (PPT helper) | ✅ | `FindAcceptorSites_PPTContribution_VerifiedByHelper` |

**Summary:** 17/17 tests, 0 missing, 0 weak, 0 duplicate.

### 5.3 Consolidation Status

Consolidation complete. Old acceptor tests removed from shared file.

| File | Role | Test Count |
|------|------|------------|
| `SpliceSitePredictor_AcceptorSite_Tests.cs` | Canonical SPLICE-ACCEPTOR-001 tests | 17 |
| `SpliceSitePredictorTests.cs` | Shared tests for other splice Test Units | 19 |

---

## 6. Assumption Register

**Total assumptions:** 0

All previously documented assumptions have been resolved:

| # | Former Assumption | Resolution | Evidence |
|---|-------------------|------------|----------|
| 1 | PWM weights approximate | **Verified** — values at key positions (-3: C=0.70, -2: A=1.00, -1: G=1.00, 0: G=0.50) match Shapiro & Senapathy (1987) ranges. Upstream pyrimidine enrichment (C+U=0.80) matches documented 70-80% | Shapiro & Senapathy (1987) |
| 2 | Normalization formula implementation-specific | **Design decision** — heuristic linear normalization mapping composite PWM + PPT scores to [0, 1]. Behavioural properties (range, monotonicity, relative ordering) verified by tests M5, M6, S5 | Tests M5, M6, S5 |
| 3 | U12 acceptor fixed 0.6 score | **Resolved** — replaced with YCCAC consensus scoring per Hall & Padgett (1994) and Jackson (1991). Scores based on match to Y-C-C-A-C pattern + PPT quality | Hall & Padgett (1994), Minor spliceosome Wikipedia |

---

## 7. Bug Fixes

| Date | Bug | Fix | Impact |
|------|-----|-----|--------|
| 2026-03-16 | AcceptorPwm alignment off-by-2 in `ScoreAcceptorSite` — PWM offsets are splice-site-relative but were applied to A-position | Changed `int pos = position + offset` → `int pos = position + 2 + offset` | Canonical site with strong PPT: score 0.27 → 0.84. Gene structure prediction now detects proper introns. |

---

## 8. Open Questions / Decisions

None.
