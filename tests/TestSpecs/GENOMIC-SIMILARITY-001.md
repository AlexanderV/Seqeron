# Test Specification: GENOMIC-SIMILARITY-001

**Test Unit ID:** GENOMIC-SIMILARITY-001
**Area:** Analysis
**Algorithm:** Sequence Similarity (k-mer Jaccard index)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Jaccard P. (1901). Étude comparative de la distribution florale… Bull. Soc. Vaudoise Sci. Nat. 37(142):547–579 (primary, via Wikipedia) | 1 (via rank-4) | https://en.wikipedia.org/wiki/Jaccard_index | 2026-06-14 |
| 2 | Ondov et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology 17:132 | 1 | https://doi.org/10.1186/s13059-016-0997-x · https://pmc.ncbi.nlm.nih.gov/articles/PMC4915045/ | 2026-06-14 |
| 3 | Mash documentation — Distance Estimation | 3 | https://github.com/marbl/Mash/blob/master/doc/sphinx/distances.rst | 2026-06-14 |

### 1.2 Key Evidence Points

1. Jaccard index: `J(A,B) = |A∩B| / |A∪B|`, the intersection size over the union size — Source 1.
2. Range `0 ≤ J(A,B) ≤ 1`; J=1 for identical sets, J=0 for disjoint sets — Source 1.
3. Applied to k-mer sets, J is "the fraction of shared k-mers out of all distinct k-mers in A and B" — Source 2.
4. k-mers compared as **sets** (distinct hashes); repeated k-mers counted once — Source 2/3.
5. J is undefined when the union is empty ("not well-defined when μ(A∪B)=0") — Source 1.

### 1.3 Documented Corner Cases

- Identical sets → J=1 (max); disjoint sets → J=0 (min) — Source 1.
- Empty union (both sets empty) → Jaccard undefined; implementation returns 0 (ASSUMPTION-1).
- Distinct-k-mer set semantics: within-sequence repeats counted once — Source 2.

### 1.4 Known Failure Modes / Pitfalls

1. Treating k-mers as a multiset/bag instead of a set would change results — the cited definition uses distinct k-mers — Source 2.
2. Empty-union division by zero — undefined per Source 1; guarded by implementation (returns 0).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateSimilarity(DnaSequence, DnaSequence, int kmerSize=5)` | GenomicAnalyzer | Canonical | k-mer Jaccard index ×100 |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Result lies in [0, 100] (Jaccard ∈ [0,1] scaled ×100) | Yes | Source 1 (0 ≤ J ≤ 1) |
| INV-2 | Identical non-empty sequences → exactly 100.0 (J=1) | Yes | Source 1 |
| INV-3 | Disjoint k-mer sets → 0.0 (J=0) | Yes | Source 1 |
| INV-4 | Symmetric: `f(a,b,k) = f(b,a,k)` | Yes | Source 1 (set symmetry) |
| INV-5 | k-mers compared as sets (within-sequence repeats counted once) | Yes | Source 2 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Partial overlap, exact fraction | ACGTACGT vs ACGTACGA, k=3; A={ACG,CGT,GTA,TAC}, B=A∪{CGA}; 4/5 | 80.0 | Src 1,2 |
| M2 | Identical sequences | ACGTACGT vs itself, k=3; J=1 | 100.0 | Src 1 |
| M3 | Disjoint k-mer sets | AAAAA vs CCCCC, k=3; {AAA} vs {CCC}; 0/2 | 0.0 | Src 1 |
| M4 | Non-integer fraction | ACGT vs ACGA, k=3; {ACG,CGT} vs {ACG,CGA}; 1/3 | 100.0/3 (33.333…) | Src 1 |
| M5 | Distinct-set semantics | AAAAAA vs AAAA, k=3; both = {AAA}; repeats counted once; 1/1 | 100.0 | Src 2 |
| M6 | Symmetry | f(ACGTACGT, ACGTACGA, 3) == f(ACGTACGA, ACGTACGT, 3) | equal (both 80.0) | Src 1 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Both empty (empty union) | "" vs "", k=3; union empty | 0.0 | ASM-1 convention |
| S2 | Both shorter than k | AC vs GT, k=3; both k-mer sets empty | 0.0 | ASM-1 convention |
| S3 | Null sequence1 | null first arg | ArgumentNullException | failure mode |
| S4 | Null sequence2 | null second arg | ArgumentNullException | failure mode |
| S5 | Invalid kmerSize | kmerSize=0 | ArgumentOutOfRangeException | failure mode |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Range invariant | varied inputs all in [0,100] | 0 ≤ r ≤ 100 | INV-1 |
| C2 | One side empty | ACGTAC vs "", k=3 | 0.0 | one set empty, other non-empty |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for existing similarity tests. `GenomicAnalyzerTests.cs` and the `GenomicAnalyzer_*` files contain no `CalculateSimilarity` test. No prior coverage exists for this unit.

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
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_CalculateSimilarity_Tests.cs` — all cases for this unit.
- **Remove:** none (no pre-existing similarity tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| GenomicAnalyzer_CalculateSimilarity_Tests.cs | Canonical | 13 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented evidence-based test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented (100/3 exact) | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented (set semantics) | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented (symmetry) | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented | ✅ Done |
| 11 | S5 | ❌ Missing | Implemented | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented (range invariant) | ✅ Done |
| 13 | C2 | ❌ Missing | Implemented (one side empty) | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact 80.0 |
| M2 | ✅ Covered | 100.0 |
| M3 | ✅ Covered | 0.0 |
| M4 | ✅ Covered | 100.0/3 |
| M5 | ✅ Covered | set semantics |
| M6 | ✅ Covered | symmetry |
| S1 | ✅ Covered | empty/empty → 0 |
| S2 | ✅ Covered | short/short → 0 |
| S3 | ✅ Covered | ArgumentNullException |
| S4 | ✅ Covered | ArgumentNullException |
| S5 | ✅ Covered | ArgumentOutOfRangeException |
| C1 | ✅ Covered | range [0,100] |
| C2 | ✅ Covered | one side empty → 0 |

Total in-scope cases: 13. ✅ count: 13.

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| ASM-1 | Empty-union returns 0.0 (Jaccard undefined for empty union) | S1, S2, C2 |
| ASM-2 | ×100 percentage scaling (presentation only) | all numeric cases |
| ASM-3 | Default k=5 (resolution only; tests pass k explicitly) | none directly (k passed explicitly) |

---

## 7. Open Questions / Decisions

1. Suffix tree NOT used: the metric is a set intersection/union of distinct k-mers, not an occurrence-search problem. A `HashSet<string>` gives O(n+m) construction and O(1) membership — optimal here; a suffix tree adds construction cost without a matching query. Decision recorded in algorithm doc §5.2.
2. Empty-union value (0 vs 1) is undefined by the cited sources; implementation contract = 0 (ASM-1).
