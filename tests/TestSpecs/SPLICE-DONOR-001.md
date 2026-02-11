# Test Specification: SPLICE-DONOR-001

**Test Unit ID:** SPLICE-DONOR-001
**Area:** Splicing
**Algorithm:** Donor (5') Splice Site Detection
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

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
4. Position -1 (G) ~80% conserved, -2 (A) ~60%, -3 (A/C) ~35% each — S1
5. GC-AG introns are valid (~0.5-1% of U2-type) — S2
6. U12-type AT-AC introns use AT at donor (~0.3%) — S2
7. PWM log-odds scoring is standard for splice site prediction — S3
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
| `ScoreDonorSite(context)` | SpliceSitePredictor | **Internal** | Private; tested indirectly via FindDonorSites |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned site has Type == Donor (or U12Donor if non-canonical) | Yes | Algorithm design |
| INV-2 | Score is in [0, 1] range | Yes | Normalization formula |
| INV-3 | Confidence is in [0, 1] range | Yes | CalculateConfidence clamps |
| INV-4 | Empty/null input yields no results | Yes | Guard clause |
| INV-5 | Sequence < 6 chars yields no results | Yes | Guard clause |
| INV-6 | Perfect consensus (CAGGUAAGU) scores ≥ any weaker context | Yes | PWM properties — S1 |
| INV-7 | GC donor score < equivalent GT donor score | Yes | 0.7 penalty — **ASSUMPTION** |
| INV-8 | Motif string is non-empty for every returned site | Yes | GetMotifContext logic |
| INV-9 | All sites are at valid positions within sequence bounds | Yes | Loop bounds |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | CanonicalGT_ConsensusMotif | CAGGUAAGU (perfect consensus) produces a donor site with Score > 0 | ≥1 site with Type=Donor, Score > 0 | S1: GU invariant at +0/+1 |
| M2 | NoGU_ReturnsEmpty | AAAAACCCCC has no GU → empty result | Empty collection | Trivially correct |
| M3 | EmptyInput_ReturnsEmpty | Empty string → empty result | Empty collection | Trivially correct |
| M4 | ShortSequence_ReturnsEmpty | Sequence < 6 chars → empty result | Empty collection | Implementation guard |
| M5 | StrongVsWeak_ScoreOrdering | CAGGUAAGU scores higher than UUUGUAAUU | strong.Score > weak.Score | S1, S3: PWM log-odds |
| M6 | GC_NonCanonical_Detected | CAGGCAAGU with includeNonCanonical=true finds GC donor | ≥1 site found | S2: GC-AG introns valid |
| M7 | GC_NotDetected_WhenCanonicalOnly | CAGGCAAGU with includeNonCanonical=false → no GT present → empty | Empty or no canonical donor | S2 |
| M8 | DNA_T_Equivalence | CAGGTAAGT (T) produces same site count as CAGGUAAGU (U) | Same number of sites | Implementation T→U conversion |
| M9 | LowercaseHandling | cagguaagu lowercase → finds donor site | ≥1 site | Implementation ToUpperInvariant |
| M10 | MultipleSites_AllDetected | Sequence with 2 separate GT contexts → finds ≥2 | ≥2 sites | Algorithm scans full sequence |

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
| C1 | U12_AT_Donor | AU donor with includeNonCanonical → detected as U12Donor | U12Donor type site | Minor spliceosome |
| C2 | BetaGlobin_RealSequence | HBB intron 1 donor context → detected | At least one donor site | Published gene |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Existing tests in `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictorTests.cs`
- Donor-specific tests: lines 12-76 (6 tests) + lines 430-458 (2 input handling tests)

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 — Canonical GT detection | ⚠ Weak | `FindDonorSites_CanonicalGT_FindsSite`: asserts ≥1 and any Donor type, no score check |
| M2 — No GU returns empty | ✅ Covered | `FindDonorSites_NoGT_ReturnsEmpty` |
| M3 — Empty input | ❌ Missing | Not tested |
| M4 — Short sequence | ✅ Covered | `FindDonorSites_ShortSequence_ReturnsEmpty` |
| M5 — Score ordering | ❌ Missing | No comparative scoring test |
| M6 — GC non-canonical | ⚠ Weak | `FindDonorSites_NonCanonicalGC_WhenEnabled`: only checks count, no type/score |
| M7 — GC not detected canonical-only | ❌ Missing | Not tested as separate case |
| M8 — DNA T equivalence | ⚠ Weak | `FindDonorSites_HandlesDNA_T`: only checks Is.Not.Null |
| M9 — Lowercase handling | ⚠ Weak | `FindDonorSites_HandlesLowercase`: only checks Is.Not.Null |
| M10 — Multiple sites | ⚠ Weak | `FindDonorSites_MultipleGT_FindsAll`: asserts ≥2 but no position check |
| S1 — Score range | ❌ Missing | Only in integration test (not donor-specific) |
| S2 — Confidence range | ❌ Missing | Only in integration test |
| S3 — Motif non-empty | ⚠ Weak | `FindDonorSites_ReturnsMotifContext`: conditional check |
| S5 — Null input | ❌ Missing | Not tested for FindDonorSites |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs` — all SPLICE-DONOR-001 tests
- **Remove from shared file:** Donor Site Tests region, donor-related input handling tests
- **Keep in shared file:** Acceptor, branch point, intron, gene structure, alternative splicing, MaxEnt, coding region, integration tests (for future Test Units)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SpliceSitePredictor_DonorSite_Tests.cs` | Canonical for SPLICE-DONOR-001 | ~17 |
| `SpliceSitePredictorTests.cs` | Remaining tests for SPLICE-ACCEPTOR/PREDICT | Reduced |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | PWM values are approximations of Shapiro & Senapathy statistics | INV-6, M5 |
| A2 | Score normalization formula is implementation-specific | INV-2, S1 |
| A3 | GC donor 0.7 penalty is an implementation heuristic | INV-7, M6 |

---

## 7. Open Questions / Decisions

None. Evidence is sufficient to define correct behavior for all MUST tests.
