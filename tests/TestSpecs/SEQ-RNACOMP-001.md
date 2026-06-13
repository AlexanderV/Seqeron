# Test Specification: SEQ-RNACOMP-001

**Test Unit ID:** SEQ-RNACOMP-001
**Area:** Composition
**Algorithm:** RNA-specific Complement (per-base, IUPAC-complete)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Cornish-Bowden A. (1985). NC-IUB nomenclature for incompletely specified bases (recommendations 1984). *Nucleic Acids Research* 13(9):3021–3030. | 2 | https://doi.org/10.1093/nar/13.9.3021 | 2026-06-13 |
| 2 | Biopython `Bio/Data/IUPACData.py` — `ambiguous_rna_complement` | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/IUPACData.py | 2026-06-13 |
| 3 | Biopython `Bio/Seq.py` — `complement_rna` / `_rna_complement_table` | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py | 2026-06-13 |
| 4 | Biopython 1.79 `Bio.Seq` docs — `complement_rna`, `reverse_complement_rna` examples | 3 | https://biopython.org/docs/1.79/api/Bio.Seq.html | 2026-06-13 |
| 5 | Bioinformatics.org SMS — IUPAC codes table | 5 | https://www.bioinformatics.org/sms/iupac.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. RNA complement table (verbatim, Biopython `ambiguous_rna_complement`): A→U, C→G, G→C, U→A, M→K, R→Y, W→W, S→S, Y→R, K→M, V→B, H→D, D→H, B→V, X→X, N→N — Source 2.
2. T is treated as U in the RNA complement (T→A), built as `ambiguous_rna_complement["T"]=...["U"]` — Sources 3, 4 ("Any T in the sequence is treated as a U").
3. Forward worked example: `complement_rna("ACG")` → `"UGC"` — Source 4.
4. Full-alphabet worked example: `reverse_complement_rna("ACGTUacgtuXYZxyz")` → `"zrxZRXaacguAACGU"`; un-reversed forward complement = `"UGCAAugcaaXRZxrz"` — Source 4.
5. Code↔complement is identical to DNA except U replaces T in the emitted alphabet; self-complementary W,S,N (and X); reciprocal A↔U, C↔G, R↔Y, M↔K, D↔H, B↔V — Sources 2, 5.

### 1.3 Documented Corner Cases

- T in RNA context → A (absorbed into U) — Sources 3, 4.
- X → X (explicit map); characters absent from the table (Z, gaps `.`/`-`, digits) pass through unchanged — Sources 2, 5.

### 1.4 Known Failure Modes / Pitfalls

1. Using the DNA complement (A→T) instead of the RNA one (A→U) — distinguished by Source 4 example `complement_rna("ACG")="UGC"`.
2. Failing to treat T as U in an RNA sequence (would yield T→A as in DNA, which is coincidentally the same letter, but the conceptual rule is T≡U) — Sources 3, 4.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetRnaComplementBase(char)` | `SequenceExtensions` | Canonical | Static, O(1) per base; IUPAC-complete RNA complement. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Canonical RNA pairing: A↔U and C↔G are reciprocal | Yes | Sources 2, 4, 5 |
| INV-2 | Output alphabet is RNA: no recognized base maps to `T`; A maps to `U` (not `T`) | Yes | Sources 3, 4 |
| INV-3 | T is treated as U → `GetRnaComplementBase('T') == 'A'` | Yes | Sources 3, 4 |
| INV-4 | Self-complementary codes: W→W, S→S, N→N (and X→X) | Yes | Sources 2, 5 |
| INV-5 | Reciprocal ambiguity pairs: R↔Y, M↔K, D↔H, B↔V | Yes | Sources 2, 5 |
| INV-6 | Recognized bases return uppercase regardless of input case (repo convention) | Yes | ASSUMPTION (SEQ-COMP-001 MUST-02); §6 |
| INV-7 | Non-IUPAC characters (gap, digit, Z) pass through unchanged | Yes | Sources 2, 5 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Standard pairing uppercase | A,U,C,G complements | A→U, U→A, C→G, G→C | Sources 2,4,5 |
| M2 | T treated as U | `GetRnaComplementBase('T')` | 'A' | Sources 3,4 |
| M3 | IUPAC ambiguity codes (uppercase) | R,Y,S,W,K,M,B,D,H,V,N | Y,R,S,W,M,K,V,H,D,B,N | Sources 2,5 |
| M4 | Lowercase recognized input → uppercase | a,u,c,g,t,r,y,n | U,A,G,C,A,Y,R,N | Source 2 mapping + repo convention (SEQ-COMP-001 MUST-02) |
| M5 | Full-alphabet worked example | each char of `"ACGTUacgtuXYZxyz"` | Biopython forward complement `"UGCAAugcaaXRZxrz"`; under repo convention recognized bases (incl. Y/y→R) uppercase and non-IUPAC chars (X, Z, x, z) pass through verbatim → `"UGCAAUGCAAXRZxRz"` | Source 4 |
| M6 | RNA-specific vs DNA | `GetRnaComplementBase('A')` differs from `GetComplementBase('A')` | 'U' vs 'T' | Source 4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Non-IUPAC pass-through | gap `-`, `.`, digit `5`, letter `Z`, space | returned unchanged | Sources 2,5 |
| S2 | Involution on RNA alphabet | complement(complement(x)) for x in {A,U,C,G,R,Y,S,W,K,M,B,D,H,V,N} | returns x | INV-1,4,5 (Sources 2,5) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | No recognized base emits T | every recognized base/code | result != 'T' | INV-2; emphasizes RNA alphabet |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `GetRnaComplementBase`. No existing test references this method (sibling `SequenceExtensions_Complement_Tests.cs` tests only the DNA `GetComplementBase`). New unit → all planned cases start ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceExtensions_GetRnaComplementBase_Tests.cs` — all cases for this unit.
- **Remove:** none (no pre-existing tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceExtensions_GetRnaComplementBase_Tests.cs` | Canonical fixture for SEQ-RNACOMP-001 | 9 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented (T→A) | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented all 11 ambiguity codes | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented lowercase→uppercase | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented full-string worked example | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented RNA-vs-DNA distinction | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented pass-through | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented involution property | ✅ Done |
| 9 | C1 | ❌ Missing | Implemented no-T-emitted | ✅ Done |

**Total items:** 9
**✅ Done:** 9 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `Standard..._ReturnsRnaComplement` |
| M2 | ✅ Covered | `Thymine_ReturnsAdenine` |
| M3 | ✅ Covered | `IupacAmbiguityCodes_..._RnaComplements` |
| M4 | ✅ Covered | `LowercaseRecognized_ReturnsUppercase` |
| M5 | ✅ Covered | `BiopythonAlphabetExample_..._MatchesPerBase` |
| M6 | ✅ Covered | `Adenine_DiffersFromDnaComplement` |
| S1 | ✅ Covered | `NonIupacCharacters_PassThroughUnchanged` |
| S2 | ✅ Covered | `Involution_AllRnaBasesAndCodes` |
| C1 | ✅ Covered | `RecognizedBases_NeverEmitThymine` |

Total in-scope cases: 9. ✅ count: 9.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Recognized bases return uppercase regardless of input case (repo convention, mirrors SEQ-COMP-001 DNA sibling; Biopython preserves case). Identity of complement is unaffected. | INV-6, M4, M5 |

---

## 7. Open Questions / Decisions

1. Decision: Mirror the established repository case convention (uppercase recognized bases) rather than Biopython case-preservation; recorded as Assumption #1 and in algorithm doc §5.4. No other open questions — complement identities are fully source-backed.
