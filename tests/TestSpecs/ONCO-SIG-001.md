# Test Specification: ONCO-SIG-001

**Test Unit ID:** ONCO-SIG-001
**Area:** Oncology
**Algorithm:** SBS-96 Single-Base-Substitution Trinucleotide Context Catalog (pyrimidine-strand folding)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Alexandrov et al. (2013), Nature 500:415-421 | 1 | https://www.nature.com/articles/nature12477 | 2026-06-14 |
| 2 | COSMIC SBS96 (Sanger) | 2/5 | https://cancer.sanger.ac.uk/signatures/sbs/sbs96/ | 2026-06-14 |
| 3 | Bergstrom et al. (2019), SigProfilerMatrixGenerator, BMC Genomics 20:685 | 1/3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6717374/ | 2026-06-14 |
| 4 | Complementarity (molecular biology), Wikipedia | 4 | https://en.wikipedia.org/wiki/Complementarity_(molecular_biology) | 2026-06-14 |

### 1.2 Key Evidence Points

1. The six substitution subtypes are C>A, C>G, C>T, T>A, T>C, T>G, each referred to by the pyrimidine of the
   mutated Watson-Crick base pair — COSMIC SBS96; Alexandrov (2013).
2. 96 mutation types = 6 substitutions × 4 5'-bases × 4 3'-bases; the mutated base is centred in the
   trinucleotide — COSMIC SBS96; SigProfiler (Bergstrom 2019).
3. Mutations whose reference (mutated) base is a purine (A or G) are folded to the pyrimidine strand by taking
   the reverse complement of the trinucleotide context and the substitution — SigProfiler (Bergstrom 2019):
   "using the purine base ... will require taking the reverse complement sequence".
4. Complement map A↔T, C↔G — Watson-Crick pairing (source 4).

### 1.3 Documented Corner Cases

- Purine reference base must be reverse-complemented before counting (SigProfiler).
- Only single-base substitutions are SBS-96 events; indels/DBS/MBS are excluded (SigProfiler/COSMIC).
- Non-ACGT context base, ref==alt: not classifiable (derived from the definition).

### 1.4 Known Failure Modes / Pitfalls

1. Counting a purine-reference mutation without folding produces an invalid (non-pyrimidine) channel —
   SigProfiler (Bergstrom 2019).
2. Including non-substitution variants inflates SBS channels — COSMIC/SigProfiler classification scope.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifySbsContext(fivePrime, refBase, altBase, threePrime)` | OncologyAnalyzer | Canonical | Folds one SBS to its 96-channel pyrimidine label |
| `Build96ContextCatalog(variants)` | OncologyAnalyzer | Canonical | Tallies SBS variants into the 96 channels |
| `EnumerateSbs96Channels()` | OncologyAnalyzer | Canonical | Returns the 96 canonical channel labels |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every classified channel uses a pyrimidine (C or T) reference base | Yes | COSMIC SBS96; SigProfiler |
| INV-2 | The channel space has exactly 96 distinct labels (6×4×4) | Yes | COSMIC SBS96; Alexandrov 2013 |
| INV-3 | Σ catalog counts = number of classifiable SBS variants (partition) | Yes | Definition (catalogue is a partition) |
| INV-4 | Folding is reverse-complement: purine-ref input maps to the revcomp pyrimidine channel | Yes | SigProfiler (Bergstrom 2019) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Pyrimidine C>A unchanged | A,C,A,A | "A[C>A]A" | COSMIC SBS96 |
| M2 | Pyrimidine C>T unchanged | T,C,T,G | "T[C>T]G" | COSMIC SBS96 |
| M3 | Pyrimidine T>C unchanged | G,T,C,A | "G[T>C]A" | COSMIC SBS96 |
| M4 | Purine G>T fold | T,G,T,A | "T[C>A]A" | SigProfiler revcomp rule + worked example |
| M5 | Purine A>G fold | C,A,G,T | "A[T>C]G" | SigProfiler revcomp rule + worked example |
| M6 | Purine G>C fold | G,G,C,C | "G[C>G]C" | SigProfiler revcomp rule + worked example |
| M7 | Purine A>T fold | A,A,T,A | "T[T>A]T" | SigProfiler revcomp rule + worked example |
| M8 | Enumerate 96 channels | EnumerateSbs96Channels() | 96 distinct labels, all pyrimidine-ref, covering all 6 subs × 16 contexts | COSMIC (6×4×4); INV-2 |
| M9 | Catalog counts | tally a known multiset incl. a folded purine variant | per-channel counts match hand tally; Σ=#SBS | INV-3; folding rule |
| M10 | Catalog folds purine into same channel as its pyrimidine form | T[G>T]A and T[C>A]A co-counted | both increment "T[C>A]A" | INV-4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null variants | Build96ContextCatalog(null) | ArgumentNullException | input validation |
| S2 | Empty variants | Build96ContextCatalog([]) | all 96 channels present with count 0; Σ=0 | partition of empty set |
| S3 | Non-ACGT context base | ClassifySbsContext('N','C','A','A') | ArgumentException | no defined context |
| S4 | ref == alt | ClassifySbsContext('A','C','C','A') | ArgumentException | not a substitution |
| S5 | Multi-char / invalid base | ClassifySbsContext('A','X','A','A') | ArgumentException | not A/C/G/T |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Lower-case input | ClassifySbsContext('a','c','a','a') | "A[C>A]A" | case-insensitive robustness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing SBS / trinucleotide / 96-context code or tests in `OncologyAnalyzer.cs` or the test project
  (grep for trinucleotide|sbs|96|signature returned only unrelated z-score `1.96` literals). This is a
  brand-new unit; every planned case starts as ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10 | ❌ Missing | new unit |
| S1–S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifySbsContext_Tests.cs` —
  all cases for the three methods, `#region` per method.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_ClassifySbsContext_Tests.cs | canonical | 16 |

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
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | S1 | ❌ Missing | implemented | ✅ Done |
| 12 | S2 | ❌ Missing | implemented | ✅ Done |
| 13 | S3 | ❌ Missing | implemented | ✅ Done |
| 14 | S4 | ❌ Missing | implemented | ✅ Done |
| 15 | S5 | ❌ Missing | implemented | ✅ Done |
| 16 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact channel label |
| M2 | ✅ Covered | exact channel label |
| M3 | ✅ Covered | exact channel label |
| M4 | ✅ Covered | exact folded label |
| M5 | ✅ Covered | exact folded label |
| M6 | ✅ Covered | exact folded label |
| M7 | ✅ Covered | exact folded label |
| M8 | ✅ Covered | 96 distinct pyrimidine labels |
| M9 | ✅ Covered | hand-tallied counts + sum |
| M10 | ✅ Covered | purine + pyrimidine co-counted |
| S1 | ✅ Covered | ArgumentNullException |
| S2 | ✅ Covered | 96 zero-count channels |
| S3 | ✅ Covered | ArgumentException |
| S4 | ✅ Covered | ArgumentException |
| S5 | ✅ Covered | ArgumentException |
| C1 | ✅ Covered | case-insensitive |

**In-scope cases:** 16 | **✅:** 16

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Channel label rendered as `5'[REF>ALT]3'` (display form; the partition is identical to any rendering) | channel keys |

---

## 7. Open Questions / Decisions

1. Channel vector ordering (the order of the 96 labels in a vector) is a presentation detail and not
   correctness-affecting for per-variant classification; the catalog is keyed by explicit label, so order does
   not change which variant falls in which channel. `EnumerateSbs96Channels()` returns them in
   substitution-major (C>A,C>G,C>T,T>A,T>C,T>G) then 5'(A,C,G,T) then 3'(A,C,G,T) order for determinism.
