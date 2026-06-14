# Test Specification: SEQ-ATSKEW-001

**Test Unit ID:** SEQ-ATSKEW-001
**Area:** Composition
**Algorithm:** AT Skew — (A − T) / (A + T)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. Mol Biol Evol 13(5):660–665. | 1 | https://doi.org/10.1093/oxfordjournals.molbev.a025626 | 2026-06-14 |
| 2 | Charneski et al. (2011). Atypical AT Skew in Firmicute Genomes. PLoS Genet 7(9):e1002283. | 1 | https://doi.org/10.1371/journal.pgen.1002283 | 2026-06-14 |
| 3 | Wikipedia "GC skew" (citing Lobry 1996). | 4 | https://en.wikipedia.org/wiki/GC_skew | 2026-06-14 |
| 4 | Biopython `Bio.SeqUtils.GC_skew` source. | 3 | https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/__init__.py | 2026-06-14 |

### 1.2 Key Evidence Points

1. AT skew = **(A − T) / (A + T)** — Charneski et al. (2011) abstract; corroborated by Wikipedia citing Lobry (1996).
2. Skew value range is [−1, +1]; AT skew = +1 ⇔ T = 0, AT skew = −1 ⇔ A = 0 — Wikipedia/Lobry (1996).
3. Base counting is case-insensitive (`count("G")+count("g")` analog) — Biopython `GC_skew` source.
4. Zero denominator (A + T = 0) ⇒ skew = 0.0 (ZeroDivisionError caught) — Biopython `GC_skew` source.
5. Non-A/T symbols (G, C, N, gaps, ambiguity) are ignored — Biopython docstring "does NOT look at any ambiguous nucleotides".

### 1.3 Documented Corner Cases

- A + T = 0 (no A and no T, e.g. all G/C) ⇒ result 0.0 (Biopython).
- Non-ACGT characters contribute to neither numerator nor denominator (Biopython).
- Pure-A ⇒ +1.0; pure-T ⇒ −1.0 (Wikipedia/Lobry bounds).

### 1.4 Known Failure Modes / Pitfalls

1. Counting non-A/T bases in the denominator would shrink the skew magnitude — incorrect; only A and T count (Biopython).
2. Returning NaN/throwing on A + T = 0 instead of 0.0 — incorrect per reference implementation (Biopython).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateAtSkew(string)` | GcSkewCalculator | **Canonical** | Core (A−T)/(A+T) over raw string; uppercases input. |
| `CalculateAtSkew(DnaSequence)` | GcSkewCalculator | **Delegate** | Forwards to the same core on the normalized `DnaSequence.Sequence`; smoke-tested for equivalence + null. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Result ∈ [−1, +1] for any input. | Yes | Wikipedia/Lobry (1996) range. |
| INV-2 | A = T (and A + T > 0) ⇒ result = 0. | Yes | Formula (A−T)/(A+T) = 0. |
| INV-3 | A + T = 0 ⇒ result = 0 (no exception, no NaN). | Yes | Biopython ZeroDivisionError → 0.0. |
| INV-4 | Result is case-insensitive (lowercase ≡ uppercase). | Yes | Biopython case-insensitive counting; repo ToUpperInvariant. |
| INV-5 | Symbols other than A/T do not change the value. | Yes | Biopython "does NOT look at ambiguous nucleotides". |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Pure A | `"AAAA"` (A=4,T=0) | 1.0 | Bounds: T=0 ⇒ +1 (Wikipedia/Lobry) |
| M2 | Pure T | `"TTTT"` (A=0,T=4) | −1.0 | Bounds: A=0 ⇒ −1 (Wikipedia/Lobry) |
| M3 | Balanced | `"ATAT"` (A=2,T=2) | 0.0 | (A−T)/(A+T)=0 (Charneski 2011) |
| M4 | Asymmetric positive | `"AAAT"` (A=3,T=1) | 0.5 | (3−1)/4 (Charneski 2011) |
| M5 | Asymmetric negative | `"ATTT"` (A=1,T=3) | −0.5 | (1−3)/4 (Charneski 2011) |
| M6 | No A/T | `"GGCC"` (A=0,T=0) | 0.0 | Zero denominator ⇒ 0.0 (Biopython) |
| M7 | G/C ignored | `"AAATGGGCCC"` (A=3,T=1) | 0.5 | Only A/T counted (Biopython) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Lowercase | `"aaat"` | 0.5 | Case-insensitive (Biopython) |
| S2 | Null string | `(string)null` | 0.0 | Documented validation (returns 0) |
| S3 | Empty string | `""` | 0.0 | Documented validation (returns 0) |
| S4 | Null DnaSequence | `(DnaSequence)null` | ArgumentNullException | Documented validation |
| S5 | Range bound | sequence with mixed A/T magnitude stays within [−1,1] | within [−1,1] | INV-1 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Overload equivalence | DnaSequence overload == string overload | equal | Delegate proof |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `GcSkewCalculator.CalculateAtSkew`. Implementation present in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs` (region "AT Skew Calculation"); no tests referenced it.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M7 | ❌ Missing | No prior tests. |
| S1–S5 | ❌ Missing | No prior tests. |
| C1 | ❌ Missing | No prior tests. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GcSkewCalculator_CalculateAtSkew_Tests.cs` — all M/S/C cases for both overloads.
- **Remove:** none (no prior tests existed).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GcSkewCalculator_CalculateAtSkew_Tests.cs` | Canonical (this unit) | 13 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented | ✅ Done |
| 12 | S5 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `CalculateAtSkew_PureAdenine_ReturnsPlusOne` |
| M2 | ✅ Covered | `CalculateAtSkew_PureThymine_ReturnsMinusOne` |
| M3 | ✅ Covered | `CalculateAtSkew_BalancedAt_ReturnsZero` |
| M4 | ✅ Covered | `CalculateAtSkew_ExcessAdenine_ReturnsExactFraction` |
| M5 | ✅ Covered | `CalculateAtSkew_ExcessThymine_ReturnsNegativeFraction` |
| M6 | ✅ Covered | `CalculateAtSkew_NoAdenineOrThymine_ReturnsZero` |
| M7 | ✅ Covered | `CalculateAtSkew_IgnoresGcAndOtherSymbols` |
| S1 | ✅ Covered | `CalculateAtSkew_LowercaseInput_MatchesUppercase` |
| S2 | ✅ Covered | `CalculateAtSkew_NullString_ReturnsZero` |
| S3 | ✅ Covered | `CalculateAtSkew_EmptyString_ReturnsZero` |
| S4 | ✅ Covered | `CalculateAtSkew_NullDnaSequence_Throws` |
| S5 | ✅ Covered | `CalculateAtSkew_AnyInput_StaysWithinBounds` |
| C1 | ✅ Covered | `CalculateAtSkew_DnaSequenceOverload_MatchesStringOverload` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Lowercase + non-A/T symbol handling for the AT-skew analog taken from Biopython `GC_skew` (no AT-skew-specific source line). The formula itself is fully sourced. | S1, M7, INV-4, INV-5 |

---

## 7. Open Questions / Decisions

1. None. Formula and range are fully sourced (Charneski 2011, Wikipedia/Lobry 1996); symbol-handling convention is sourced by analogy to Biopython `GC_skew` and matches the existing implementation.
