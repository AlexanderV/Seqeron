# Test Specification: SEQ-COMPLEX-COMPRESS-001

**Test Unit ID:** SEQ-COMPLEX-COMPRESS-001
**Area:** Complexity
**Algorithm:** Lempel–Ziv complexity (compression-based sequence complexity)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Lempel & Ziv (1976), On the Complexity of Finite Sequences, IEEE TIT 22(1):75–81 | 1 | https://doi.org/10.1109/TIT.1976.1055501 | 2026-06-14 |
| 2 | Wikipedia, Lempel–Ziv complexity (cites #1) | 4 | https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity | 2026-06-14 |
| 3 | Naereen/Lempel-Ziv_Complexity (Python reference, MIT) | 3 | https://github.com/Naereen/Lempel-Ziv_Complexity/blob/master/src/lempel_ziv_complexity.py | 2026-06-14 |
| 4 | entropy/antropy `lziv_complexity` (cites #1, #5) | 3 | https://raphaelvallat.com/entropy/build/html/generated/entropy.lziv_complexity.html | 2026-06-14 |
| 5 | Zhang et al. (2009), Normalized LZ complexity, J Math Chem 46(4):1203–1212 | 1 | https://doi.org/10.1007/s10910-008-9512-2 | 2026-06-14 |

### 1.2 Key Evidence Points

1. LZ complexity = number of distinct substrings (components) produced parsing the sequence left-to-right; a new component starts where the running substring is no longer a previously-encountered word — source #1/#2/#3.
2. Reference parser (set-based exhaustive history): grow the running substring while it is already in the seen-set; otherwise add it and restart — source #3.
3. Worked exact values: `1001111011000010`→8, `1010101010101010`→7, `1001111011000010000010`→9, `100111101100001000001010`→10 — source #3 doctests.
4. Normalization: `LZn = c / (n / log_b n)` with `b` = alphabet size (distinct symbols) — source #4 (citing #5).
5. Asymptotic upper bound `b(n) = n/log_α(n)`; normalized value → 1 for random sequences — source #6 (WebSearch synthesis of primary-citing papers).

### 1.3 Documented Corner Cases

- Empty sequence → complexity 0 (no components) — traced reference parser.
- Homopolymer `"0"×16` → components `0/00/000/0000/00000` → c=5 (`c = ⌊(√(8n+1)−1)/2⌋`) — traced reference parser.
- Single distinct symbol (b<2) → `log_b n` undefined → normalization returns the raw count — source #4 rule.

### 1.4 Known Failure Modes / Pitfalls

1. Wrong parsing convention: the older `entropy` doc parses `1001111011000010` into 6 components; the LZ76 exhaustive-history convention (this unit) yields 8 — sources #3 vs #4. We follow the exhaustive-history convention (sources #1/#2/#3).
2. Trailing-partial-component counting differs by ±1 between the Wikipedia pseudocode and the set-based reference; we adopt the set-based contract (source #3) — see ASSUMPTION A1.
3. `log` base must be alphabet size, not 2 nor e, for normalization — source #4.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateLempelZivComplexity(DnaSequence)` | SequenceComplexity | **Canonical** | Raw LZ76 component count |
| `CalculateLempelZivComplexity(string)` | SequenceComplexity | **Canonical** | Raw LZ76 component count (string) |
| `CalculateNormalizedLempelZivComplexity(DnaSequence)` | SequenceComplexity | **Canonical** | `c/(n/log_b n)` |
| `CalculateNormalizedLempelZivComplexity(string)` | SequenceComplexity | **Canonical** | normalized (string) |
| `EstimateCompressionRatio(DnaSequence)` | SequenceComplexity | **Delegate** | returns normalized LZ complexity |
| `EstimateCompressionRatio(string)` | SequenceComplexity | **Delegate** | returns normalized LZ complexity |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Empty/null-string input → raw complexity 0 | Yes | traced reference parser (source #3) |
| INV-2 | Raw complexity ≥ 1 for any non-empty sequence | Yes | each first symbol is a component (source #3) |
| INV-3 | Raw complexity ≤ n for length-n input | Yes | each component is ≥1 char (source #1/#3) |
| INV-4 | A homopolymer has strictly lower complexity than a string of all-distinct symbols of the same length | Yes | productivity buildup (source #1/#2) |
| INV-5 | `EstimateCompressionRatio` equals `CalculateNormalizedLempelZivComplexity` for the same input | Yes | delegation (design) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Doctest 1 | `CalculateLempelZivComplexity("1001111011000010")` | 8 | source #3 doctest |
| M2 | Doctest 2 | `CalculateLempelZivComplexity("1010101010101010")` | 7 | source #3 doctest |
| M3 | Doctest 3 | `CalculateLempelZivComplexity("1001111011000010000010")` | 9 | source #3 doctest |
| M4 | Doctest 4 | `CalculateLempelZivComplexity("100111101100001000001010")` | 10 | source #3 doctest |
| M5 | Homopolymer | `CalculateLempelZivComplexity("0000000000000000")` | 5 | traced parser (source #3 rule) |
| M6 | All-distinct | `CalculateLempelZivComplexity("ACGT")` | 4 | parsing rule (source #3) |
| M7 | Normalized | `CalculateNormalizedLempelZivComplexity("1001111011000010")` | 2.0 (8/(16/log₂16)) | source #4 formula + derivation |
| M8 | b<2 fallback | `CalculateNormalizedLempelZivComplexity("0000000000000000")` | 5.0 (raw count) | source #4 undefined-log rule |
| M9 | Delegation | `EstimateCompressionRatio("1001111011000010")` equals normalized (2.0) | 2.0 | INV-5 (design) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty string | `CalculateLempelZivComplexity("")` | 0 | INV-1 |
| S2 | Null DnaSequence | `CalculateLempelZivComplexity((DnaSequence)null)` | ArgumentNullException | sibling convention |
| S3 | Single base | `CalculateLempelZivComplexity("A")` | 1 | INV-2 |
| S4 | DnaSequence overload parity | `"ACGT"` via DnaSequence | 4 | overload consistency |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-4 property | homopolymer < all-distinct of same length | true | invariant |
| C2 | DNA normalization (b=4) | normalized uses log base 4 | matches `c/(n/log₄ n)` | Zhang 2009 application |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexityTests.cs` — contained 4 pre-existing `EstimateCompressionRatio` tests asserting the OLD non-source-backed heuristic (exact 14/27, 5/112, range [0,1]) plus an empty/null guard. The heuristic-output and [0,1]-range tests are invalid under the corrected Lempel–Ziv implementation.
- `tests/Seqeron/Seqeron.Mcp.Sequence.Tests/ComplexityCompressionRatioTests.cs` — MCP binding test asserting `CompressionRatio <= 1` (old heuristic range) and a repetitive-vs-diverse ordering using a single-symbol low example.
- No canonical `SequenceComplexity_EstimateCompressionRatio_Tests.cs` existed prior to this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| M8 | ❌ Missing | new unit |
| M9 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_EstimateCompressionRatio_Tests.cs` — all M/S/C cases for this unit.
- **Remove:** the heuristic-output `EstimateCompressionRatio_HighComplexity_ReturnsExact` (14/27), `EstimateCompressionRatio_LowComplexity_ReturnsExact` (5/112), and `EstimateCompressionRatio_RangeIsZeroToOne` tests from `SequenceComplexityTests.cs` (they asserted the replaced heuristic / an invalid [0,1] bound). Kept: `EstimateCompressionRatio_EmptySequence_ReturnsZero`, `EstimateCompressionRatio_NullSequence_ThrowsException`.
- **Fix:** `Seqeron.Mcp.Sequence.Tests/ComplexityCompressionRatioTests.cs` — drop the invalid `<= 1` bound (normalized LZ may exceed 1) and make the repetitive-vs-diverse ordering compare two length-40 four-letter sequences so the log_b(n) factor is held constant.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceComplexity_EstimateCompressionRatio_Tests.cs` | canonical | 15 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | M9 | ❌ Missing | implemented | ✅ Done |
| 10 | S1 | ❌ Missing | implemented | ✅ Done |
| 11 | S2 | ❌ Missing | implemented | ✅ Done |
| 12 | S3 | ❌ Missing | implemented | ✅ Done |
| 13 | S4 | ❌ Missing | implemented | ✅ Done |
| 14 | C1 | ❌ Missing | implemented | ✅ Done |
| 15 | C2 | ❌ Missing | implemented | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact doctest value 8 |
| M2 | ✅ Covered | exact doctest value 7 |
| M3 | ✅ Covered | exact doctest value 9 |
| M4 | ✅ Covered | exact doctest value 10 |
| M5 | ✅ Covered | homopolymer 5 |
| M6 | ✅ Covered | all-distinct 4 |
| M7 | ✅ Covered | normalized 2.0 |
| M8 | ✅ Covered | b<2 fallback 5.0 |
| M9 | ✅ Covered | delegation 2.0 |
| S1 | ✅ Covered | empty → 0 |
| S2 | ✅ Covered | null → ArgumentNullException |
| S3 | ✅ Covered | single base → 1 |
| S4 | ✅ Covered | DnaSequence parity |
| C1 | ✅ Covered | INV-4 property |
| C2 | ✅ Covered | DNA b=4 normalization |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Trailing-partial-component convention follows the set-based reference (source #3), not the +1 of the Wikipedia pseudocode | parsing contract, M1–M6 |
| A2 | Normalization log base = number of distinct symbols actually present (b); b<2 returns raw count | M7, M8, C2 |

---

## 7. Open Questions / Decisions

1. Decision: `EstimateCompressionRatio` (the registry-canonical name) is retained as a thin delegate returning the normalized LZ complexity, replacing the prior non-source-backed heuristic. Raw and normalized LZ are exposed as new canonical methods.
2. Decision: exhaustive-history (LZ76) parsing convention adopted over the alternative `entropy` convention because it matches the primary description and has reproducible worked values.
