# Test Specification: RNA-INVERT-001

**Test Unit ID:** RNA-INVERT-001
**Area:** RnaStructure
**Algorithm:** RNA Inverted Repeats (potential stem regions)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Alamro et al. 2021, IUPACpal (BMC Bioinformatics) | 1 | https://doi.org/10.1186/s12859-021-03983-2 (full text https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/) | 2026-06-14 |
| 2 | Ussery, Wassenaar & Borini 2008 (via Wikipedia "Inverted repeat") | 4 | https://en.wikipedia.org/wiki/Inverted_repeat | 2026-06-14 |
| 3 | EMBOSS einverted manual (Rice et al. 2000) | 3 | https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted | 2026-06-14 |

### 1.2 Key Evidence Points

1. An inverted repeat (IR) is a sequence followed downstream by its reverse complement — formally `WGW̄ᴿ`, W = left arm, G = gap/loop (|G| ≥ 0), W̄ᴿ = right arm = reverse complement of W. — IUPACpal (PMC7866733); Wikipedia/Ussery 2008.
2. The right arm equals the **reverse complement** (not just complement, not parallel) of the left arm; pairing is antiparallel. — IUPACpal; Wikipedia example `TTACGnnnnnnCGTAA` (CGTAA = revcomp(TTACG)).
3. A perfect IR is the k = 0 (Hamming-distance 0) case of the gapped-IR-within-k-mismatches model. — IUPACpal.
4. Inverted repeats are equivalent to potential stem-loops; found by comparing the sequence to its reverse complement; arms are complementary, the intervening region is the loop. — EMBOSS einverted.
5. Watson-Crick complement basis A⟷U, C⟷G (A⟷T for DNA). — IUPACpal.

### 1.3 Documented Corner Cases

- Loop length |G| ≥ 0; |G| = 0 yields a palindrome (IUPACpal). The default RNA call requires a loop (minSpacing ≥ 3).
- Arms shorter than the minimum length are not reported; arms extend to their maximal complementary run (IUPACpal; einverted maximal local alignment).

### 1.4 Known Failure Modes / Pitfalls

1. Parallel-direct repeats (left arm == right arm read 5'→3') are NOT inverted repeats — the right arm must be the *reverse* complement. — IUPACpal / Wikipedia ("reverse complement").
2. Treating G-U wobble as complementary would over-report; this repository uses strict Watson-Crick + IUPAC complement only. — `GetRnaComplementBase`.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindInvertedRepeats(string sequence, int minLength = 4, int minSpacing = 3, int maxSpacing = 100)` | RnaSecondaryStructure | **Canonical** | Returns `(int Start1, int End1, int Start2, int End2, int Length)` tuples: left arm `[Start1..End1]`, right arm `[Start2..End2]`, arm `Length`. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | For every reported repeat, the right arm equals the reverse complement of the left arm under strict Watson-Crick/IUPAC: for k in [0,Length), `complement(seq[Start2+Length-1-k]) == seq[Start1+k]`. | Yes | IUPACpal WGW̄ᴿ; `GetRnaComplementBase` |
| INV-2 | `End1 - Start1 + 1 == Length` and `End2 - Start2 + 1 == Length`; both arms have equal length. | Yes | IUPACpal (W and W̄ᴿ same length) |
| INV-3 | The loop length `Start2 - End1 - 1` lies in `[minSpacing, maxSpacing]`. | Yes | IUPACpal max gap; einverted loop |
| INV-4 | `Length >= minLength`; arms shorter than minLength are not reported. | Yes | IUPACpal minimum length |
| INV-5 | `Start1 < End1 < Start2 < End2` (arms are disjoint, left precedes right). | Yes | IUPACpal WGW̄ᴿ ordering |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Known IR, exact positions | `UUACGAAAAAACGUAA`, defaults | exactly one repeat `(0,4,11,15,5)` (left UUACG, right CGUAA = revcomp) | Wikipedia/Ussery `TTACGnnnnnnCGTAA`; IUPACpal |
| M2 | Palindromic IR, 3-nt loop | `GGCCAAAGGCC`, defaults | exactly one repeat `(0,3,7,10,4)` | IUPACpal WGW̄ᴿ |
| M3 | Antiparallel required (parallel rejected) | `AAGGAAAAAGG` (arms AAGG==AAGG parallel) | empty (AAGG is not revcomp of AAGG) | IUPACpal/Wikipedia "reverse complement" |
| M4 | No inverted repeat | `AAAAAAAAAAAA` (no complementary arm) | empty | definition (no WGW̄ᴿ) |
| M5 | Reverse-complement invariant holds | result of M1/M2 satisfies INV-1 antiparallel pairing for every k | all pairs complementary | IUPACpal; INV-1 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Arm shorter than minLength | `GGCCAAAGGCC`, minLength = 5 | empty (arm only 4) | INV-4 boundary |
| S2 | Loop below minSpacing | `GGCCAGGCC` (loop length 1), defaults (minSpacing 3) | empty | INV-3 lower bound |
| S3 | Loop above maxSpacing | `UUACGAAAAAACGUAA`, maxSpacing = 5 (loop is 6) | empty | INV-3 upper bound |
| S4 | Null / empty input | `null`, `""` | empty (no throw) | sibling-method convention (yield break) |
| S5 | Too-short input | `GGCC` (< minLength*2+minSpacing) | empty | length guard |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Maximal-arm extension | `GGGGCAAAGCCCC` (5-nt perfect arm) | arm Length 5, not truncated to 4 | einverted maximal local alignment |
| C2 | Property: invariants on random-ish input | for several fixed sequences, every reported repeat satisfies INV-1..INV-5 | all hold | O(n²) invariant property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- An implementation existed: `RnaSecondaryStructure.FindInvertedRepeats` (Analysis project). No dedicated test file (`RnaSecondaryStructure_FindInvertedRepeats_Tests.cs`) existed. The existing implementation located the right arm using `j = i + minLength + spacing` (gap measured from `minLength`, not the matched arm length), so it could not correctly locate arms longer than `minLength` — Present-but-nonconforming (corrected in Phase 5).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | no prior test |
| M2 | ❌ Missing | no prior test |
| M3 | ❌ Missing | no prior test |
| M4 | ❌ Missing | no prior test |
| M5 | ❌ Missing | no prior test |
| S1 | ❌ Missing | no prior test |
| S2 | ❌ Missing | no prior test |
| S3 | ❌ Missing | no prior test |
| S4 | ❌ Missing | no prior test |
| S5 | ❌ Missing | no prior test |
| C1 | ❌ Missing | no prior test |
| C2 | ❌ Missing | no prior test |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_FindInvertedRepeats_Tests.cs` — all M/S/C cases.
- **Remove:** nothing (no prior test file for this unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `RnaSecondaryStructure_FindInvertedRepeats_Tests.cs` | Canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-position test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented exact-position test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented parallel-rejection test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented no-IR test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented revcomp invariant test | ✅ Done |
| 6 | S1 | ❌ Missing | Implemented minLength boundary test | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented minSpacing test | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented maxSpacing test | ✅ Done |
| 9 | S4 | ❌ Missing | Implemented null/empty test | ✅ Done |
| 10 | S5 | ❌ Missing | Implemented too-short test | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented maximal-arm test | ✅ Done |
| 12 | C2 | ❌ Missing | Implemented invariant property test | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact tuple `(0,4,11,15,5)` |
| M2 | ✅ | exact tuple `(0,3,7,10,4)` |
| M3 | ✅ | empty (parallel rejected) |
| M4 | ✅ | empty |
| M5 | ✅ | INV-1 verified per-position |
| S1 | ✅ | empty at minLength 5 |
| S2 | ✅ | empty (loop < minSpacing) |
| S3 | ✅ | empty (loop > maxSpacing) |
| S4 | ✅ | empty, no throw |
| S5 | ✅ | empty (length guard) |
| C1 | ✅ | arm Length 5 |
| C2 | ✅ | INV-1..INV-5 hold |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Perfect (k=0), ungapped arms only (strict WC/IUPAC complement); scored mismatch/gap variant out of scope | M1–M5, C1 |
| 2 | minSpacing/maxSpacing bound the loop length |G| | M1, M2, S2, S3 |
| 3 | Maximal arm + non-overlapping greedy reporting | M1, C1 |

---

## 7. Open Questions / Decisions

1. Decision: the existing implementation's right-arm offset bug (`j = i + minLength + spacing`) is corrected to extend from the matched arm length; behavior now conforms to the WGW̄ᴿ definition. No remaining open questions.
