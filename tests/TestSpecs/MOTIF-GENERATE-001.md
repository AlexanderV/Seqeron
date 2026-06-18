# Test Specification: MOTIF-GENERATE-001

**Test Unit ID:** MOTIF-GENERATE-001
**Area:** Matching
**Algorithm:** IUPAC-Degenerate Consensus Generation (`MotifFinder.GenerateConsensus`)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Cornish-Bowden / NC-IUB (1985). Nomenclature for incompletely specified bases. NAR 13(9):3021. | 2 | https://doi.org/10.1093/nar/13.9.3021 | 2026-06-14 |
| 2 | UCSC Genome Browser — IUPAC ambiguity codes | 5 | https://genome.ucsc.edu/goldenPath/help/iupac.html | 2026-06-14 |
| 3 | Wikipedia — Nucleic acid notation (Table 1, cites NC-IUB 1984) | 4 | https://en.wikipedia.org/wiki/Nucleic_acid_notation | 2026-06-14 |
| 4 | DECIPHER `ConsensusSequence` (Bioconductor) | 3 | https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. IUPAC set→symbol mapping is bijective: {A,G}=R, {C,T}=Y, {C,G}=S, {A,T}=W, {G,T}=K, {A,C}=M, {C,G,T}=B, {A,G,T}=D, {A,C,T}=H, {A,C,G}=V, {A,C,G,T}=N — source [1][2][3].
2. A degenerate consensus combines, at each column, the bases that pass a frequency threshold into the IUPAC symbol for that base set; minority bases below the threshold are removed — source [4].
3. When ≥2 bases pass the inclusion rule the column emits a degeneracy code, not a single base — source [4].
4. N denotes all four bases; a single missing base yields a three-base not-X code (B/D/H/V), not N — source [1][2].

### 1.3 Documented Corner Cases

- Frequency threshold governs base inclusion before IUPAC encoding (DECIPHER [4]).
- Equal-abundance bases → degeneracy code (DECIPHER [4]).
- This implementation's threshold is `count > total × 0.25` (strict `>`), a documented design constant; the 25 % value is implementation-specific (DECIPHER's own default differs), see §6.

### 1.4 Known Failure Modes / Pitfalls

1. Treating a base at exactly the threshold as included — boundary is strict `>`, so exactly-25 % bases are excluded (this implementation).
2. Emitting N for four-equal columns — under strict `>` 25 % no base passes, so the fallback most-frequent base is emitted, not N (implementation contract; see §6).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GenerateConsensus(IEnumerable<string>)` | `MotifFinder` | **Canonical** | IUPAC-degenerate consensus; threshold = count > n×0.25 |
| `GetIupacCode(...)` | `MotifFinder` (private) | **Internal** | set→symbol mapping; tested indirectly via `GenerateConsensus` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Output length equals the length of the first input sequence | Yes | Column-wise per-position construction [4] |
| INV-2 | A unanimous column (single base) yields that standard base (A/C/G/T) | Yes | Singleton set → standard base [1] |
| INV-3 | Each output character is one of the 15 IUPAC symbols {A,C,G,T,R,Y,S,W,K,M,B,D,H,V,N} | Yes | NC-IUB symbol alphabet [1][2] |
| INV-4 | The symbol emitted for a passing base set is exactly the NC-IUB symbol for that set | Yes | Bijective mapping [1][2][3] |
| INV-5 | A base with count ≤ n×0.25 is excluded from the code (strict `>` boundary) | Yes | Implementation design constant; threshold family [4] |
| INV-6 | Empty input collection → empty string | Yes | Guard contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | TwoBase_AG_R | column {A,G} both passing | `"R"` | NC-IUB [1] |
| M2 | TwoBase_CT_Y | column {C,T} | `"Y"` | NC-IUB [1] |
| M3 | TwoBase_CG_S | column {C,G} | `"S"` | NC-IUB [1] |
| M4 | TwoBase_AT_W | column {A,T} | `"W"` | NC-IUB [1] |
| M5 | TwoBase_GT_K | column {G,T} | `"K"` | NC-IUB [1] |
| M6 | TwoBase_AC_M | column {A,C} | `"M"` | NC-IUB [1] |
| M7 | ThreeBase_CGT_B | column {C,G,T} all passing | `"B"` | NC-IUB [1] |
| M8 | ThreeBase_AGT_D | column {A,G,T} | `"D"` | NC-IUB [1] |
| M9 | ThreeBase_ACT_H | column {A,C,T} | `"H"` | NC-IUB [1] |
| M10 | ThreeBase_ACG_V | column {A,C,G} | `"V"` | NC-IUB [1] |
| M11 | Unanimous_ReturnsInput | identical sequences | input string verbatim | INV-2 [1] |
| M12 | MultiColumn_MixedCodes | `["ATGC","GTGC"]` col0={A,G}→R, rest unanimous | `"RTGC"` | NC-IUB [1] |
| M13 | ThresholdBoundary_Exactly25Excluded | `["AAAA","AAGT","AACT","AATT"]` col3 T(2)>1.0, others ≤1.0 | col3 = `'T'` | INV-5 (design constant) |
| M14 | MinorityBelowThreshold_Dropped | `["AAGGC"]→` split as A,A,G,G,C col; C(1)≤1.25 dropped | `"R"` | DECIPHER threshold [4] |
| M15 | NoBasePasses_FallbackMostFrequent | `["A","C","G","T"]` none >1.0 → most-frequent, tie→A | `"A"` | implementation contract §6 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | CaseInsensitive_LowerUpper | lowercase input | same as upper-cased | input normalisation |
| S2 | Empty_ReturnsEmpty | empty collection | `""` | INV-6 |
| S3 | OutputLength_MatchesFirst | length invariant | length == first seq length | INV-1 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null_Throws | null collection | `ArgumentNullException` | guard |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/.../MotifFinderTests.cs` (region "Consensus Sequence Tests"): 5 tests for `GenerateConsensus` — several use permissive assertions.
- `tests/.../MutationKillerTests.cs` (region "MotifFinder — Consensus generation survivors"): 3 mutation-killer tests for `GenerateConsensus`.
- No canonical `MotifFinder_GenerateConsensus_Tests.cs` existed.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| MotifFinderTests.GenerateConsensus_IdenticalSequences_ReturnsSame | 🔁 Duplicate | replaced by M11 |
| MotifFinderTests.GenerateConsensus_MixedBases_ReturnsIupac | ⚠ Weak | `.Or.EqualTo('A')/('G')` — alternative outcomes; replaced by M1/M12 |
| MotifFinderTests.GenerateConsensus_Empty_ReturnsEmpty | 🔁 Duplicate | replaced by S2 |
| MotifFinderTests.GenerateConsensus_AllDifferent_ReturnsMostCommon | ⚠ Weak | `Does.Match("^[ACGTN]+$")` — shape not value; replaced by M15 |
| MotifFinderTests.GenerateConsensus_NullSequences_ThrowsException | 🔁 Duplicate | replaced by C1 |
| MutationKillerTests.GenerateConsensus_ThresholdBoundary_ExactlyAtQuarter | 🔁 Duplicate | replaced by M13 (exact, full positions) |
| MutationKillerTests.GenerateConsensus_NoPresentBases_FallbackToMaxBy | ⚠ Weak | `BeOneOf('A','C','G','T')` — replaced by M15 (exact 'A') |
| MutationKillerTests.GenerateConsensus_TwoBases_ReturnsAmbiguityCode | 🔁 Duplicate | replaced by M1 |
| M1–M15, S1–S3, C1 | ❌ Missing | implement in canonical file |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_GenerateConsensus_Tests.cs` — all evidence-based cases for `GenerateConsensus`.
- **Remove:** the `GenerateConsensus_*` tests in `MotifFinderTests.cs` (region "Consensus Sequence Tests" + the null test) and the "MotifFinder — Consensus generation survivors" region in `MutationKillerTests.cs` — duplicated/weakened by the canonical file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MotifFinder_GenerateConsensus_Tests.cs` | Canonical | 19 |
| `MotifFinderTests.cs` | Other MotifFinder methods (consensus tests removed) | (unchanged for other methods) |
| `MutationKillerTests.cs` | Other mutation survivors (consensus region removed) | (unchanged for other regions) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1–M10 | ❌ Missing | implemented exact set→symbol tests | ✅ Done |
| 2 | M11 | ❌ Missing | implemented unanimous test | ✅ Done |
| 3 | M12 | ❌ Missing | implemented multi-column mixed | ✅ Done |
| 4 | M13 | ❌ Missing | implemented exact-25% boundary (full positions) | ✅ Done |
| 5 | M14 | ❌ Missing | implemented minority-dropped | ✅ Done |
| 6 | M15 | ❌ Missing | implemented fallback most-frequent (exact 'A') | ✅ Done |
| 7 | S1–S3 | ❌ Missing | implemented case/empty/length | ✅ Done |
| 8 | C1 | ❌ Missing | implemented null guard | ✅ Done |
| 9 | MotifFinderTests consensus tests | 🔁/⚠ | removed | ✅ Done |
| 10 | MutationKillerTests consensus region | 🔁/⚠ | removed | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M15 | ✅ | canonical file, exact values |
| S1–S3 | ✅ | canonical file |
| C1 | ✅ | canonical file |
| Old MotifFinderTests consensus tests | ✅ | removed (consolidated) |
| Old MutationKillerTests consensus region | ✅ | removed (consolidated) |

All in-scope cases ✅. Count of ✅ = total in-scope cases.

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | 25 % strict-`>` inclusion threshold is a documented design constant (threshold-consensus family is authoritative; exact 25 % is implementation-specific) | M13, M14, M15, INV-5 |
| 2 | Fallback to single most-frequent base (alphabetical tie-break) when no base passes the threshold | M15 |
| 3 | Length from first sequence; case-insensitive; non-ACGT ignored | INV-1, S1 |

---

## 7. Open Questions / Decisions

1. The 25 % threshold is correctness-affecting but documented and named in code; the *symbol* output for any given passing base set is dictated by the authoritative NC-IUB table, which is fully source-backed. Tests pin the boundary explicitly and otherwise use unambiguous inputs so verified symbols depend only on the authoritative table. No unresolved correctness-affecting assumption blocks completion.
