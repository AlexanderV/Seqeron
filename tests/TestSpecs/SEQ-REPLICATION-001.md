# Test Specification: SEQ-REPLICATION-001

**Test Unit ID:** SEQ-REPLICATION-001
**Area:** Composition
**Algorithm:** Replication Origin Prediction (cumulative GC-skew minimum)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-10-09

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Grigoriev A (1998), Nucleic Acids Res 26(10):2286–2290 | 1 | https://doi.org/10.1093/nar/26.10.2286 | 2026-06-14 |
| 2 | Lobry JR (1996), Mol Biol Evol 13(5):660–665 | 1 | https://pubmed.ncbi.nlm.nih.gov/8676740/ | 2026-06-14 |
| 3 | Rosalind, Minimum Skew Problem (BA1F) | 3 | https://rosalind.info/problems/ba1f/ | 2026-06-14 |
| 4 | Wikipedia, GC skew (cited primaries 1,2) | 4 | https://en.wikipedia.org/wiki/GC_skew | 2026-06-14 |
| 5 | Rosalind BA1F extra dataset (go-rosalind mirror) | 3 | https://raw.githubusercontent.com/charlesreid1/go-rosalind/master/rosalind/data/minimum_skew.txt | 2026-10-09 |
| 6 | Lu & Salzberg (2020) SkewIT + `src/skewi.py` | 1 | https://doi.org/10.1371/journal.pcbi.1008439 ; https://github.com/jenniferlu717/SkewIT | 2026-10-09 |

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
| `FindMinimumSkewPositions` / `FindMaximumSkewPositions(DnaSequence\|string, bool circular = false)` | GcSkewCalculator | Canonical | All BA1F minimizers / maximizers, ascending |
| `PredictReplicationOrigin(DnaSequence\|string, bool circular)` | GcSkewCalculator | Canonical | Positions mod n (n ≡ 0) |
| `PredictReplicationOrigin(DnaSequence\|string, int windowSize)` | GcSkewCalculator | Canonical | Grigoriev windowed cumulative diagram extrema (window centres) |
| `CalculateSkewIndex(DnaSequence\|string, int windowSize = 20000)` | GcSkewCalculator | Canonical | SkewIT SkewI (`skewi.py`), null when skewi.py reports none |

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

### 4.2b Finisher additions (2026-10-09, A1-2/A1-3) — `GcSkewCalculator_ReplicationOriginExtensions_Tests`, `GcSkewCalculator_SkewIndex_Tests`

| ID | Test Case | Expected Outcome | Evidence |
|----|-----------|------------------|----------|
| A1 | BA1F sample all minimizers | `[53, 97]` | Rosalind sample output |
| A2 | BA1F extra dataset (embedded) | `[89969, 89970, 89971, 90345, 90346]`; max `[20377, 20378, 20379]`; −184/+154 | published answer + python |
| A3 | BA1F sample all maximizers | `[16, 20, 21]` | python brute force |
| A4 | small ties (`CCGGCC`, `GGGCCC`, `CGCG`, `AATT`, lower case) | brute-force sets | python |
| A5 | null/empty | null DnaSequence throws; string → `[0]` | contract |
| C4 | circular mod n (`GGGCCC` → {0}, `G` → max 0, …) | brute force over Skew_0..Skew_{n−1} | python |
| C5 | Skew_n = 0 rotation-equivariance (r = 0, 12345, 50000, 99999) | (p + r) mod n = original sets | python; Grigoriev arbitrary start |
| C6 | Skew_n ≠ 0 non-equivariance (`CCGGG` vs rotation 3) | {2} vs {0, 4} | python |
| C7 | circular=false ≡ original; BA1F circular unchanged | equality | — |
| C8 | circular null/empty | throws / zero prediction | contract |
| W1 | windowed prediction (synthetic w1000; extra w1000, w5000; BA1F w10) | Biopython cumsum values (Evidence table) | Biopython 1.88 + numpy |
| W2 | window 1 on `GGGCCC` | origin 5 (0), terminus 2 (+3) | definition (no Skew_0) |
| W3 | windowed guards | w < 1 throws; short → zero | contract |
| K1 | SkewI small sequences k4 | skewi.py values (Evidence) | SkewIT code run |
| K2 | SkewI none (12 windows / no G-C) | null | SkewIT code run |
| K3 | SkewI BA1F extra k1000 / k20 | 0.203158581311549 / 0.034430033253851994 | SkewIT code run |
| K4 | SkewI ideal genome | 1.0 (cap) | SkewIT code run |
| K5 | defaults/guards/case | k=20000; 5 windows → null; lower case = upper | SkewIT code + documented deviation |
| M-W | MCP `predict_replication_origin` `windowSize` (F29) | BA1F sample w10 → 45 / −0.3095238095238095, 15 / 0.5; w25 → 37 / −0.15384615384615385, 62 / 0.11888111888111885; position arrays empty; `circular`+`windowSize` → ArgumentException; w/k < 1 → AOORE; "ACG" w10 → zero, not significant | Biopython 1.88 `numpy.cumsum(GC_skew(seq, w)[:n//w])` |
| M-K | MCP `predict_replication_origin` `skewIndexWindow` (F29) | k4: 52-mer 0.23076923076923078, 58-mer 0.6206896551724138, BA1F sample 0.16, 12 windows → null; omitted → null; independent of `windowSize` | SkewIT `skewi.py -k 4 --min-len 0` (run) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null DnaSequence | `PredictReplicationOrigin((DnaSequence)null!)` | ArgumentNullException | input validation |
| C2 | Null/empty string | `PredictReplicationOrigin((string)null!)`, `""` | zero prediction, IsSignificant=false | documented handling |
| C3 | Single base | `G` → origin/terminus | origin=0 (skew 0), terminus=1 (skew +1) | boundary |
| R1 | Synthetic genome vs python reference | `ACGTC`×6000+`AGGTC`×10000+`ACGTC`×4000 | origin 29997 (−6000), terminus 79998 (+4001) | python BA1F re-implementation |
| R2 | Rotated circular genome | same, start rotated by 50000 | origin 79997 (−4000), terminus 29998 (+6001); map back to 29997/79998 | python reference; Grigoriev "arbitrary start" |
| R3 | Agreement with windowed Grigoriev diagram | same, `CalculateCumulativeGcSkew(g,1000)` | min −10 @ centre 29500, max 20/3 @ 79500; prediction within 1 window | `numpy.cumsum(Bio.SeqUtils.GC_skew(g,1000))` |
| R4 | Reuse of canonical cumulative skew | BA1F sample, `CalculateCumulativeGcSkew(s,1)` | min/max equal OriginSkew/TerminusSkew; terminus 16 (+2) | Biopython `GC_skew(s,1)` cumsum |

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

1. Decision: BA1F asks for ALL minimizing positions; the repository API returns a single position, so the deterministic tie-break is "first (smallest) extreme index". Documented in the algorithm doc and Evidence; tested by M4. (2026-10-09: all positions now available via `FindMinimumSkewPositions` / `FindMaximumSkewPositions`, A1–A5.)
2. Decision: the windowed-cumulative `PredictReplicationOrigin(windowSize)` overload is replaced by the canonical per-nucleotide method; legacy tests asserting the windowed model are removed as nonconforming.
