# Test Specification: SEQ-COMPLEX-COMPRESS-001

**Test Unit ID:** SEQ-COMPLEX-COMPRESS-001
**Area:** Complexity
**Algorithm:** Lempel–Ziv complexity (compression-based sequence complexity)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-28 (review-2026-09: LZ78 parse replaced by true LZ76)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Lempel & Ziv (1976), On the Complexity of Finite Sequences, IEEE TIT 22(1):75–81 | 1 | https://doi.org/10.1109/TIT.1976.1055501 | 2026-06-14 |
| 2 | Wikipedia, Lempel–Ziv complexity (cites #1) | 4 | https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity | 2026-06-14 |
| 3 | Kaspar & Schuster (1987), Easily calculable measure for the complexity of spatiotemporal patterns, Phys Rev A 36:842 (LZ76 scan algorithm) | 1 | https://doi.org/10.1103/PhysRevA.36.842 | 2026-09-28 (cited; algorithm read from antropy source) |
| 4 | antropy 0.2.2 `lziv_complexity` / `_lz_complexity` (cites #1, #5; PyPI wheel source opened) | 3 | https://pypi.org/project/antropy/ | 2026-09-28 |
| 6 | Estévez-Rams et al. (2013) arXiv:1311.0546 — exhaustive-history definition + example (WebSearch snippet) | 2 | https://arxiv.org/abs/1311.0546 | 2026-09-28 |
| 7 | Naereen/Lempel-Ziv_Complexity 0.2.2 — **counter-reference**: implements LZ78 incremental parsing, not LZ76 | 3 | https://pypi.org/project/lempel_ziv_complexity/ | 2026-09-28 |
| 5 | Zhang et al. (2009), Normalized LZ complexity, J Math Chem 46(4):1203–1212 | 1 | https://doi.org/10.1007/s10910-008-9512-2 | 2026-06-14 |

### 1.2 Key Evidence Points

1. LZ76 complexity c(S) = number of components of the exhaustive history: a component starting at p is extended while it is still a substring of the text preceding its last symbol (copy start < p, overlap allowed); the first non-reproducible extension closes it; a reproducible remainder at the end is the last component — source #1/#2/#6.
2. Reference algorithm: Kaspar–Schuster scan (source #3) = antropy `_lz_complexity` (source #4) = Wikipedia pseudocode (source #2).
3. Worked exact values: `0001101001000101`→6 (0·001·10·100·1000·101, source #1 via #6 snippet); `010011101101100`→6 (source #6); `1001111011000010`→6 (1/0/01/1110/1100/0010) and normalized 1.5, `HELLO WORLD! ×4`→11 / 0.38596001132145313, `A..Z`→26 / 1.0 — antropy doctests (source #4).
4. Normalization: `LZn = c / (n / log_b n)` with `b` = alphabet size (distinct symbols) — source #4 (citing #5).
5. Asymptotic upper bound `b(n) = n/log_α(n)`; normalized value → 1 for random sequences — source #6 (WebSearch synthesis of primary-citing papers).

### 1.3 Documented Corner Cases

- Empty sequence → complexity 0 (no components) — traced reference parser.
- Homopolymer `"0"×n` (n ≥ 2) → components `0 / 0…0` (self-overlapping copy) → c=2 — definition (#1) + antropy (#4).
- Single distinct symbol (b<2) → antropy clamps the log base to 2 (`base = 2 if base < 2 else base`) and returns `c/(n/log_2 n)`. For `"0"×16` this is `2/(16/log_2 16) = 0.5`.
- n = 1 → normalized formula undefined (`log_b 1 = 0`; antropy raises ZeroDivisionError); this library returns the raw count 1 (documented convention).

### 1.4 Known Failure Modes / Pitfalls

1. **LZ78 vs LZ76 confusion (the defect fixed 2026-09):** the "set of seen phrases" parse (Naereen, source #7) is Ziv–Lempel 1978 incremental parsing: `1001111011000010`→8 (1/0/01/11/10/110/00/010), `"0"×16`→5. The LZ76 exhaustive history gives 6 and 2. The pre-2026-09 TestSpec had these reversed.
2. The trailing reproducible remainder IS counted (Kaspar–Schuster `if len ≠ 1 then c += 1`; e.g. `010011101101100` ends with factor `101100`).
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
| M1 | antropy doctest | `CalculateLempelZivComplexity("1001111011000010")` | 6 | source #4 doctest; #2 |
| M2 | LZ 1976 example | `CalculateLempelZivComplexity("0001101001000101")` | 6 | source #1 (via #6 snippet) |
| M3 | Estévez-Rams example | `CalculateLempelZivComplexity("010011101101100")` | 6 | source #6 |
| M4 | Period 2 | `CalculateLempelZivComplexity("1010101010101010")` | 3 (1/0/10…) | antropy (#4) |
| M5 | Homopolymer | `CalculateLempelZivComplexity("0000000000000000")` | 2 | definition #1; antropy #4 |
| M5b | Text doctests | `HELLO WORLD! ×4` → 11 / 0.38596001132145313; `A..Z` → 26 / 1.0 | exact | antropy doctests (#4) |
| M6 | All-distinct | `CalculateLempelZivComplexity("ACGT")` | 4 | definition |
| M7 | Normalized | `CalculateNormalizedLempelZivComplexity("1001111011000010")` | 1.5 (6/(16/log₂16)) | antropy doctest (#4) |
| M8 | b<2 clamp | `CalculateNormalizedLempelZivComplexity("0000000000000000")` | 0.5 (`2/(16/log₂16)`) | source #4 code: `base = 2 if base < 2 else base` |
| M9 | Delegation | `EstimateCompressionRatio("1001111011000010")` equals normalized (1.5) | 1.5 | INV-5 (design) |
| R1 | DNA reference dataset | 6 DNA strings (random seed 2026, n = 20–200; CAG×20) | raw + normalized = antropy 0.2.2 | antropy (#4) |
| R2 | Long DNA | 20 000-base LCG string | c = 2756, norm = 0.984423382950957 | antropy (#4) |

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
| C2 | DNA normalization (b=4) | `ACGT×4` = A/C/G/T/ACGTACGTACGT, c = 5 → 5/(16/log₄16) = 0.625 | 0.625 | antropy (#4) |

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
| M8 | ✅ Covered | b<2 clamp-to-2 → 1.25 |
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
| A1 | (retired 2026-09) Trailing reproducible component is counted, per LZ76 definition and Kaspar–Schuster/antropy — no longer an assumption | M1–M6 |
| A2 | Normalization log base = number of distinct symbols actually present (b); b<2 is clamped to 2 (per antropy `base = 2 if base < 2 else base`); n = 1 returns the raw count 1 (formula undefined) | M7, M8, C2 |

---

## 7. Open Questions / Decisions

1. Decision: `EstimateCompressionRatio` (the registry-canonical name) is retained as a thin delegate returning the normalized LZ complexity, replacing the prior non-source-backed heuristic. Raw and normalized LZ are exposed as new canonical methods.
2. Decision (corrected 2026-09-28): true LZ76 exhaustive history (Lempel & Ziv 1976; Kaspar–Schuster 1987; antropy), replacing the LZ78-style set parse (Naereen) that had been mislabelled as LZ76. Computed via the Longest-Previous-Factor array (Crochemore & Ilie 2008) in O(n log² n), value-identical to the Kaspar–Schuster scan.
