# Test Specification: SEQ-REPLICATION-001

**Test Unit ID:** SEQ-REPLICATION-001
**Area:** Composition
**Algorithm:** Replication Origin Prediction (cumulative GC-skew minimum)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Grigoriev A (1998), Nucleic Acids Res 26(10):2286–2290 | 1 | https://doi.org/10.1093/nar/26.10.2286 | 2026-06-14 |
| 2 | Lobry JR (1996), Mol Biol Evol 13(5):660–665 | 1 | https://pubmed.ncbi.nlm.nih.gov/8676740/ | 2026-06-14 |
| 3 | Rosalind, Minimum Skew Problem (BA1F) | 3 | https://rosalind.info/problems/ba1f/ | 2026-06-14 |
| 4 | Wikipedia, GC skew (cited primaries 1,2) | 4 | https://en.wikipedia.org/wiki/GC_skew | 2026-06-14 |

### 1.2 Key Evidence Points

1. Cumulative skew Skew_i = (#G − #C) over prefix Genome[0..i); G:+1, C:−1, A/T:0; Skew_0 = 0 — Rosalind BA1F.
2. Global minimum of the cumulative diagram = replication origin; global maximum = terminus — Grigoriev 1998; Wikipedia GC skew.
3. Positions are 0-based prefix indices i ∈ [0, |Genome|]; BA1F asks for ALL minimizers — Rosalind BA1F.
4. BA1F sample genome (len 100) → minimizing positions `53 97`, min value −4 — Rosalind BA1F (re-derived in session).
5. Skews switch sign at origin/terminus; leading strand is G-rich — Lobry 1996; Grigoriev 1998.

### 1.3 Documented Corner Cases

- Ties: multiple positions may share the extreme value (BA1F returns `53 97`); this unit returns the first (smallest) minimizing/maximizing index.
- Flat diagram (no net G/C asymmetry): origin/terminus not resolvable (amplitude 0).

### 1.4 Known Failure Modes / Pitfalls

1. Using per-window skew sums (windowed cumulative) instead of per-nucleotide cumulative skew fails to reproduce the BA1F worked example — Rosalind BA1F.
2. Off-by-one in prefix indexing (Skew_0 = 0 before base 0) shifts every reported position — Rosalind BA1F.
3. Invented significance threshold (`amplitude > count × 0.01`) has no authoritative basis — removed (Evidence Assumption 1).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictReplicationOrigin(DnaSequence)` | GcSkewCalculator | Canonical | Origin = min prefix, terminus = max prefix |
| `PredictReplicationOrigin(string)` | GcSkewCalculator | Delegate | Same core; null/empty → zero prediction; case-insensitive |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | PredictedOrigin is the first prefix index minimizing the cumulative skew | Yes | Rosalind BA1F |
| INV-2 | PredictedTerminus is the first prefix index maximizing the cumulative skew | Yes | Grigoriev 1998; Wikipedia |
| INV-3 | OriginSkew = min ≤ 0 ≤ max = TerminusSkew (Skew_0 = 0 is always in range) | Yes | Rosalind BA1F (prefix starts at 0) |
| INV-4 | Positions lie in [0, n] where n = sequence length | Yes | Rosalind BA1F |
| INV-5 | IsSignificant ⇔ max > min (non-zero amplitude) | Yes | **ASSUMPTION** (Evidence §Assumptions 1) |
| INV-6 | A and T bases do not change the diagram (only G/C count) | Yes | Rosalind BA1F definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | BA1F sample origin | Full BA1F genome → origin | PredictedOrigin = 53, OriginSkew = −4 | Rosalind BA1F sample `53 97` |
| M2 | Per-nt skew increments | `CCGGGG` → min=−2@2, max=+2@6 | Origin=2 (skew −2), Terminus=6 (skew +2) | BA1F definition (re-derived) |
| M3 | Terminus = global max | `GGGCCC` → diagram 0,1,2,3,2,1,0 | Terminus=3 (skew +3); Origin=0 (skew 0) | Grigoriev 1998 (max=terminus) |
| M4 | Tie-break first index | `GGCCGGCC` → min −0? compute; first minimizing prefix | first/smallest minimizing index reported | BA1F returns multiple minimizers |
| M5 | A/T ignored (INV-6) | `GATGCA` vs `GGC` skew unaffected by A/T | A/T leave cumulative unchanged | BA1F definition |
| M6 | OriginSkew ≤ 0 ≤ TerminusSkew (INV-3) | any sequence including BA1F | min ≤ 0 and max ≥ 0 | BA1F (Skew_0 = 0) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Flat / no asymmetry | `AAAATTTT` (no G/C) | Origin=Terminus=0, skews 0, IsSignificant=false | INV-5 |
| S2 | Significant flag | BA1F sample | IsSignificant=true (max −min > 0) | INV-5 |
| S3 | Case-insensitive | lowercase BA1F-style snippet vs upper | identical result | string overload uppercases |
| S4 | Positions within bounds | `GGGCCC` | 0 ≤ origin,terminus ≤ 6 | INV-4 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null DnaSequence | `PredictReplicationOrigin((DnaSequence)null!)` | ArgumentNullException | input validation |
| C2 | Null/empty string | `PredictReplicationOrigin((string)null!)`, `""` | zero prediction, IsSignificant=false | documented handling |
| C3 | Single base | `G` → origin/terminus | origin=0 (skew 0), terminus=1 (skew +1) | boundary |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculatorTests.cs` contained legacy tests for `PredictReplicationOrigin` (`PredictReplicationOrigin_FindsMinimum`, `_FindsMaximum`, `_ExactPositionsAndSkew`, `_NullSequence_ThrowsException`, `_TooShortSequence_ReturnsDefault`) that asserted against the **windowed** (per-window sum) implementation with a `windowSize` argument and the invented `0.01×count` significance threshold. These rubber-stamped the nonconforming behavior and do not reproduce the BA1F worked example.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 BA1F sample | ❌ Missing | not previously tested |
| M2 per-nt increments | ⚠ Weak | legacy `_ExactPositionsAndSkew` used windowed sums, wrong model |
| M3 terminus=max | ⚠ Weak | legacy `_FindsMaximum` used windowed model + windowSize |
| M4 tie-break | ❌ Missing | not tested |
| M5 A/T ignored | ❌ Missing | not tested |
| M6 origin≤0≤terminus | ❌ Missing | not tested |
| S1 flat | ⚠ Weak | legacy `_TooShortSequence` asserted windowed default-zero, wrong reason |
| S2 significant flag | ⚠ Weak | legacy used 0.01×count threshold (invented) |
| S3 case-insensitive | ❌ Missing | not tested |
| S4 bounds | ❌ Missing | not tested |
| C1 null DnaSequence | ✅ Covered | legacy `_NullSequence_ThrowsException` (signature changed: no windowSize) |
| C2 null/empty string | ❌ Missing | string overload is new |
| C3 single base | ❌ Missing | not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_PredictReplicationOrigin_Tests.cs` — all SEQ-REPLICATION-001 cases.
- **Remove:** the five legacy `PredictReplicationOrigin_*` tests from `GcSkewCalculatorTests.cs` (they test the removed `windowSize` overload / invented threshold). Other Gc-skew tests in that file remain.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GcSkewCalculator_PredictReplicationOrigin_Tests.cs` | Canonical (this unit) | 16 |
| `GcSkewCalculatorTests.cs` | Other GcSkewCalculator methods (legacy PRO tests removed) | (reduced) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented BA1F sample origin test | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote against per-nt model (`CCGGGG`) | ✅ Done |
| 3 | M3 | ⚠ Weak | Rewrote terminus=max (`GGGCCC`) | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented tie-break first index | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented A/T-ignored test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented origin≤0≤terminus | ✅ Done |
| 7 | S1 | ⚠ Weak | Rewrote flat-diagram test | ✅ Done |
| 8 | S2 | ⚠ Weak | Rewrote significant flag (no threshold) | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented case-insensitive test | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented bounds test | ✅ Done |
| 11 | C1 | ✅ Covered | Re-added null DnaSequence (new signature) | ✅ Done |
| 12 | C2 | ❌ Missing | Implemented null/empty string | ✅ Done |
| 13 | C3 | ❌ Missing | Implemented single-base test | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | BA1F sample origin = 53 |
| M2 | ✅ | per-nt increments verified |
| M3 | ✅ | terminus = global max |
| M4 | ✅ | first minimizing index |
| M5 | ✅ | A/T ignored |
| M6 | ✅ | origin ≤ 0 ≤ terminus |
| S1 | ✅ | flat diagram |
| S2 | ✅ | IsSignificant = true on BA1F |
| S3 | ✅ | case-insensitive |
| S4 | ✅ | positions in [0,n] |
| C1 | ✅ | null DnaSequence throws |
| C2 | ✅ | null/empty string → zero prediction |
| C3 | ✅ | single base |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | `IsSignificant ⇔ max > min` (no authoritative numeric significance threshold exists; invented `0.01×count` removed) | INV-5, S1, S2 |

---

## 7. Open Questions / Decisions

1. Decision: BA1F asks for ALL minimizing positions; the repository API returns a single position, so the deterministic tie-break is "first (smallest) extreme index". Documented in the algorithm doc and Evidence; tested by M4.
2. Decision: the windowed-cumulative `PredictReplicationOrigin(windowSize)` overload is replaced by the canonical per-nucleotide method; legacy tests asserting the windowed model are removed as nonconforming.
