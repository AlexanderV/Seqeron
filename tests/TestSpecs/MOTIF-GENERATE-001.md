# Test Specification: MOTIF-GENERATE-001

**Test Unit ID:** MOTIF-GENERATE-001
**Area:** Matching
**Algorithm:** IUPAC-Degenerate Consensus Generation (`MotifFinder.GenerateConsensus`)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-30 (review-2026-09 B05 F13)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Cornish-Bowden / NC-IUB (1985). Nomenclature for incompletely specified bases. NAR 13(9):3021. | 2 | https://doi.org/10.1093/nar/13.9.3021 | 2026-06-14 |
| 2 | UCSC Genome Browser — IUPAC ambiguity codes | 5 | https://genome.ucsc.edu/goldenPath/help/iupac.html | 2026-06-14 |
| 3 | Wikipedia — Nucleic acid notation (Table 1, cites NC-IUB 1984) | 4 | https://en.wikipedia.org/wiki/Nucleic_acid_notation | 2026-06-14 |
| 4 | DECIPHER `ConsensusSequence` (Bioconductor) | 3 | https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html | 2026-06-14 |
| 5 | Cavener D.R. (1987) NAR 15(4):1353 — degenerate consensus rules | 1 | PMID 3822832 (WebSearch snippet) | 2026-09-30 |
| 6 | Biopython 1.88 `Bio.motifs` `degenerate_consensus` (installed source) + Tutorial chapter_motifs.rst | 3 | raw.githubusercontent.com/biopython/biopython/master/Doc/Tutorial/chapter_motifs.rst | 2026-09-30 |

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
2. Four-equal columns — no base passes strict `>` 25 %; all four tie at the maximum → `N` (DECIPHER equal-abundance rule [4]; Biopython `degenerate_consensus` = N [6]). Before F13 (2026-09) this returned `A`.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GenerateConsensus(IEnumerable<string>)` | `MotifFinder` | **Canonical** | IUPAC-degenerate consensus; threshold = count > n×0.25 |
| `GenerateConsensus(IEnumerable<string>, double inclusionThreshold)` | `MotifFinder` | **Canonical** | configurable θ ∈ [0, 1]; θ = 0.25 bit-identical to the parameterless overload (B05 follow-up) |
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
| M13′ | (M13 col2) | four-way tie column | `'N'` → full `"AANT"` | [4][6] (F13) |
| M15 | FourEqualBases_ReturnsN | `["A","C","G","T"]` none >1.0, four-way tie | `"N"` | DECIPHER [4], Biopython [6] (F13) |
| M16 | NoBasePasses_TiedBasesEncoded | `["A","C","-","-"]` → `"M"`; `["A","-","-","-"]` → `"A"` | tied-max set | [4] (F13) |
| M17 | ColumnWithoutAcgt_ReturnsN | `["A-","AN","A-"]` | `"AN"` | NC-IUB N = any [1] (F13) |
| M18 | Cavener_EqualsBiopython | 23 alignments (tutorial WACVC, GBGTW, CV; rule branches; 12 random) via `GenerateCavenerConsensus` | Biopython values | [5][6] |

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
| C2 | NullElement_Throws | null row | `ArgumentException` | F13 (was NRE) |
| C3 | UnequalLengths_Throws | `["ACG","AC"]`, `["AC","ACG"]` | `ArgumentException` | Biopython MSA (F13; was silently truncated) |
| C4 | Cavener_InvalidInput | null / null row / unequal / gap | ANE / AE | shared `BuildCountMatrix` |

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
| 2 | (resolved F13) no-pass → IUPAC code of bases tied at the maximum; no A/C/G/T → N | M15–M17 |
| 3 | Equal-length rows required; case-insensitive; non-ACGT not counted but counted in n | INV-1, S1, C3 |

---

## 7. Open Questions / Decisions

1. The 25 % threshold is correctness-affecting but documented and named in code; the *symbol* output for any given passing base set is dictated by the authoritative NC-IUB table, which is fully source-backed. Tests pin the boundary explicitly and otherwise use unambiguous inputs so verified symbols depend only on the authoritative table. No unresolved correctness-affecting assumption blocks completion.

## 8. Configurable threshold tests (B05 follow-up, 2026-09-30)

| ID | Test | Evidence |
|----|------|----------|
| T-M1 | 20 random alignments × θ locked to an independent Python oracle of the rule (`ThresholdCases`); full run 700/700 | oracle `count > θ·n`, tie fallback, NC-IUB map |
| T-M2 | θ = 0.25 ≡ parameterless overload (unit + property C1, 500 random) | bit-identity requirement |
| T-M3 | hand-derived A,A,A,C,G: θ 0/0.1 → V, 0.2/0.6 → A; A,C at θ 1 → M | rule |
| T-S1 | guards: null, θ < 0, θ > 1, NaN, unequal rows | contract |
| T-P2 | property C2: base set at θ₂ ⊆ base set at θ₁ for θ₁ ≤ θ₂ | monotonicity of the cut |

## 9. DECIPHER `ConsensusSequence` tests (B05 audit group C, F28, 2026-10-01)

| ID | Test | Evidence |
|----|------|----------|
| T-D1 | 90 stratified random cases (DNA/RNA/AA × ambiguity × includeNonLetters × includeTerminalGaps × minInformation) = output of DECIPHER 3.9.4 R/C source built against R 4.3.3 + Biostrings 2.70.2 (`DecipherCases`); full run 10,000/10,000 | `MotifFinder_DecipherConsensus_Tests` |
| T-D2 | every example of `man/ConsensusSequence.Rd` (AAAT, majority, ties, terminal gaps, `.` as gap, non-letters, degeneracy) | manual + R build |
| T-D3 | RNA → `U`, ragged rows, lower case, empty set → "" | R build |
| T-D4 | hand-derived source order: A .6/C .4 → M; A .9/C .06/G .04 → M (t .95), A (t .9) | `makeConsensus` |
| T-D5 | guards as `ConsensusSequence.R`: threshold ∉ [0,1), minInformation ∉ (0,1], noConsensusChar ∉ alphabet, characters outside DNA_/RNA_/AA_ALPHABET, null | R argument checks, Biostrings |
| T-MCP | `generate_consensus` `inclusionThreshold` (default 0.25 = parameterless; 0.2 → `HHHH`); `generate_decipher_consensus` delegation | Mcp.Analysis.Tests |
