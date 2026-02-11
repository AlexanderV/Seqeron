# Test Specification: SPLICE-ACCEPTOR-001

**Test Unit ID:** SPLICE-ACCEPTOR-001
**Area:** Splicing
**Algorithm:** Acceptor Site Detection
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Shapiro & Senapathy (1987). Nucleic Acids Res 15(17):7155–7174 | 1 | — | 2026-02-12 |
| 2 | Burge, Tuschl & Sharp (1999). The RNA World, 2nd ed. CSHL Press | 1 | — | 2026-02-12 |
| 3 | Yeo & Burge (2004). J Comput Biol 11(2–3):377–394 | 1 | — | 2026-02-12 |
| 4 | Patel & Steitz (2003). Nat Rev Mol Cell Biol 4(12):960–970 | 1 | — | 2026-02-12 |
| 5 | Wikipedia — RNA splicing | 4 | https://en.wikipedia.org/wiki/RNA_splicing | 2026-02-12 |
| 6 | Wikipedia — Polypyrimidine tract | 4 | https://en.wikipedia.org/wiki/Polypyrimidine_tract | 2026-02-12 |

### 1.2 Key Evidence Points

1. 3' splice site consensus is (Y)nNCAG|G — AG dinucleotide is almost invariant — Shapiro & Senapathy (1987)
2. PPT is 15–20 nt pyrimidine-rich region upstream of AG — Lodish et al. (2004) via Wikipedia
3. PPT quality (continuous pyrimidines) determines splice site strength — Burge et al. (1999)
4. U12-type introns use AC instead of AG at 3' splice site — Patel & Steitz (2003)
5. Position -2 = A (100%), position -1 = G (100%), position -3 = C (~70%) — Shapiro & Senapathy (1987)
6. First exonic nucleotide (position 0) favors G (~50%) — Shapiro & Senapathy (1987)

### 1.3 Documented Corner Cases

1. **U12 AC acceptor**: Minor spliceosome uses AC instead of AG — Patel & Steitz (2003)
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
| INV-5 | Position = i + 1 (position after AG) | Yes | Implementation: `Position: i + 1` |
| INV-6 | Motif is non-empty for all returned sites | Yes | GetMotifContext extracts surrounding context |
| INV-7 | Empty/null input → empty result | Yes | Guard: `if (string.IsNullOrEmpty(sequence) \|\| sequence.Length < 20) yield break` |
| INV-8 | Higher minScore → subset of lower minScore results | Yes | Score filtering: `if (score >= minScore)` |
| INV-9 | Strong PPT score ≥ weak PPT score for same AG context | Yes | PPT pyrimidine fraction contributes to score — Burge et al. (1999) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Canonical AG detected | Sequence with strong PPT followed by CAG — should find ≥1 acceptor with Type=Acceptor | ≥1 site, Type=Acceptor | Shapiro & Senapathy (1987) — AG consensus |
| M2 | No AG returns empty | All-pyrimidine sequence (no AG) — no acceptor sites | Empty result | Trivially correct |
| M3 | Empty input returns empty | `""` and `null` inputs should return empty | Empty result | Guard behavior |
| M4 | Short sequence returns empty | Sequence < 20 chars | Empty result | Guard: length < 20 |
| M5 | Strong PPT > weak PPT | Continuous pyrimidine PPT scores higher than purine-interrupted context | Strong score > weak score | Burge et al. (1999) — PPT quality |
| M6 | Score range [0, 1] | All returned scores within [0, 1] | Scores ∈ [0, 1] | INV-1 — normalization |
| M7 | Confidence range [0, 1] | All returned confidences within [0, 1] | Confidence ∈ [0, 1] | INV-2 — CalculateConfidence |
| M8 | DNA T equivalence | DNA input (T instead of U) produces same results | Same site count and scores | Implementation T→U conversion |
| M9 | Case insensitivity | Lowercase input produces same results as uppercase | Same site count | Implementation ToUpperInvariant |
| M10 | Multiple AG sites | Sequence with two AG dinucleotides in scannable range → both detected | ≥2 sites | Scanning algorithm |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | U12 AC detected | includeNonCanonical=true with AC dinucleotide → U12Acceptor found | ≥1 U12Acceptor site | Patel & Steitz (2003) |
| S2 | U12 AC excluded | includeNonCanonical=false (default) with AC → no U12Acceptor sites | No U12Acceptor in result | Default parameter behavior |
| S3 | Position after AG | Returned position = index of G in AG + 1 (first exonic nucleotide) | Position verified | INV-5 |
| S4 | Motif non-empty | All returned sites have non-empty Motif containing AG context | Non-empty motif | INV-6 |
| S5 | Threshold filtering | minScore=0.8 returns subset of minScore=0.2 results | High-threshold count ≤ low-threshold count | INV-8 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | AG before position 15 not detected | AG at position < 15 in sequence not found | Not in results | Implementation scan starts at i=15 |
| C2 | PPT pyrimidine count correlates | Independent PPT scoring helper matches expected trend | Score reflects PPT quality | **ASSUMPTION** — PPT scoring formula |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- **File:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictorTests.cs`
- **Tests found:** 5 acceptor-specific tests in `#region Acceptor Site Tests`
- **Additional:** integration tests in same file also exercise FindAcceptorSites but belong to SPLICE-PREDICT-001

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (canonical AG) | ⚠ Weak | `FindAcceptorSites_CanonicalAG_FindsSite` — asserts Count ≥ 1 and type, no score check |
| M2 (no AG → empty) | ⚠ Weak | `FindAcceptorSites_NoAG_ReturnsEmpty` — correct but minimal |
| M3 (empty/null) | ❌ Missing | No test for null or empty input |
| M4 (short sequence) | ⚠ Weak | `FindAcceptorSites_ShortSequence_ReturnsEmpty` — correct but minimal |
| M5 (strong vs weak PPT) | ⚠ Weak | `FindAcceptorSites_StrongPPT_HighScore` — doesn't compare strong vs weak scores |
| M6 (score range) | ❌ Missing | |
| M7 (confidence range) | ❌ Missing | |
| M8 (DNA T equivalence) | ❌ Missing | |
| M9 (case insensitivity) | ❌ Missing | |
| M10 (multiple AG sites) | ❌ Missing | |
| S1 (U12 AC detected) | ⚠ Weak | `FindAcceptorSites_U12NonCanonical_WhenEnabled` — asserts not empty only |
| S2 (U12 AC excluded) | ❌ Missing | |
| S3 (position after AG) | ❌ Missing | |
| S4 (motif non-empty) | ❌ Missing | |
| S5 (threshold filtering) | ❌ Missing | |
| C1, C2 | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs` — NEW, all deep acceptor tests
- **Remove:** 5 acceptor tests from `SpliceSitePredictorTests.cs` (weak, replaced by evidence-based versions)
- **Keep:** Branch point, intron, gene structure, alternative splicing, MaxEntScore, IsWithinCodingRegion, and integration tests in shared file (they belong to SPLICE-PREDICT-001)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SpliceSitePredictor_AcceptorSite_Tests.cs` | Canonical SPLICE-ACCEPTOR-001 tests | 17 |
| `SpliceSitePredictorTests.cs` | Shared tests for other splice Test Units | ~30 (after removing 5 acceptor tests) |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | PWM weights approximate — implementation values may differ from Shapiro & Senapathy (1987) due to rounding | M5, C2 |
| 2 | Normalization formula `(score/(count+1) + 2) / 4` is implementation-specific | M6, C2 |
| 3 | U12 acceptor fixed 0.6 score is a simplification | S1 |

---

## 7. Open Questions / Decisions

None.
