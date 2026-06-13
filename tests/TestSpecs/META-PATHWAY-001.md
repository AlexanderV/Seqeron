# Test Specification: META-PATHWAY-001

**Test Unit ID:** META-PATHWAY-001
**Area:** Metagenomics
**Algorithm:** Metabolic Pathway Enrichment (Over-Representation Analysis, hypergeometric test)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Boyle et al. (2004), GO::TermFinder, *Bioinformatics* 20(18):3710–3715 | 1 / 3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3037731/ | 2026-06-13 |
| 2 | PNNL Proteomics Data Analysis in R/Bioconductor §8.2 ORA | 3 | https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. Over-representation p-value = upper tail of the hypergeometric distribution:
   `P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M, n−i) / C(N, n)` — Boyle 2004; PNNL §8.2.
2. Symbols: N = background size, M = pathway/gene-set size, n = query size, x = overlap — Boyle 2004; PNNL §8.2.
3. Reference `phyper(q = x−1, m = M, n = N−M, k = n, lower.tail = FALSE)` gives P(X ≥ x) — PNNL §8.2.
4. Worked example: N=8000, M/n = {400,100}, x=20 → P(X≥20) = 7.88×10⁻⁸ — PNNL §8.2.
5. Hypergeometric (sampling without replacement) chosen over binomial for accuracy — Boyle 2004.

### 1.3 Documented Corner Cases

- x = 0: empty upper sum ⇒ P = 1 (no over-representation possible) — PNNL §8.2.
- Degenerate population (N, M, or n = 0): P(X ≥ x) = 1 — derived from the formula.
- Infeasible partial tables (i > M or n−i > N−M): C(·,·) = 0, contribute 0 — sampling-without-replacement constraint.

### 1.4 Known Failure Modes / Pitfalls

1. Using P(X > x) instead of P(X ≥ x) (off-by-one in the tail) — PNNL §8.2 (`q = x−1` correction).
2. Computing C(·) directly overflows for large N; must use log-space (log-Gamma) — implementation note.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindPathwayEnrichment(queryGenes, pathwayDatabase, backgroundGenes?)` | MetagenomicsAnalyzer | Canonical | Caller supplies pathway→genes; computes hypergeometric ORA per pathway, sorted by p-value |
| `HypergeometricUpperTail(x, bigN, bigM, n)` | MetagenomicsAnalyzer | Canonical | Core right-tail probability P(X ≥ x); directly verifiable against worked examples |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ p-value ≤ 1 for all inputs | Yes | Probability axiom; formula (Boyle 2004) |
| INV-2 | x = 0 ⇒ p-value = 1 (empty upper sum) | Yes | PNNL §8.2 corner case |
| INV-3 | Degenerate population (N, M, or n ≤ 0) ⇒ p-value = 1 | Yes | Derived from formula |
| INV-4 | P(X ≥ x) is symmetric under swapping (M ↔ n) | Yes | Hypergeometric symmetry |
| INV-5 | Results are returned ascending by p-value | Yes | Boyle 2004 (term ranking by p-value) |
| INV-6 | PathwayEnrichment fields (Overlap, PathwaySize, QuerySize, BackgroundSize) equal the counts used | Yes | Implementation contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | PNNL worked example | `HypergeometricUpperTail(20, 8000, 400, 100)` | 7.884747e-8 (≈7.88×10⁻⁸) | PNNL §8.2 |
| M2 | Symmetric orientation | `HypergeometricUpperTail(20, 8000, 100, 400)` | 7.884747e-8 | PNNL §8.2 + symmetry |
| M3 | All query in pathway | `HypergeometricUpperTail(5, 10, 5, 5)` | 1/252 = 0.0039682539683 | Formula (exact) |
| M4 | Partial overlap | `HypergeometricUpperTail(1, 4, 2, 2)` | 5/6 = 0.8333333333 | Formula (exact) |
| M5 | x = 0 corner | `HypergeometricUpperTail(0, 10, 5, 5)` | 1.0 | PNNL §8.2 corner case |
| M6 | At least one | `HypergeometricUpperTail(1, 10, 5, 5)` | 251/252 = 0.9960317460 | Formula (exact) |
| M7 | Degenerate population | `HypergeometricUpperTail(2, 0, 5, 5)` and N/M/n=0 variants | 1.0 | INV-3 / formula |
| M8 | FindPathwayEnrichment end-to-end | query {g1,g2,g3}, pathway P1={g1,g2,g3,g4,g5}, explicit background of 10 genes | Overlap=3, PathwaySize=5, QuerySize=3, BackgroundSize=10, p = HypergeometricUpperTail(3,10,5,3) = 1/12 | Formula (exact) |
| M9 | Sorting + multi-pathway | two pathways, one enriched, one not | results ascending by p-value (INV-5) | Boyle 2004 |
| M10 | Field correctness | M8 result fields | match counts (INV-6) | Implementation contract |
| M11 | Null query | `FindPathwayEnrichment(null, db)` | ArgumentNullException | Contract |
| M12 | Null pathway DB | `FindPathwayEnrichment(query, null)` | ArgumentNullException | Contract |
| M13 | Empty pathway DB | `FindPathwayEnrichment(query, {})` | empty result list | Contract |
| M14 | No overlap pathway | query disjoint from pathway | Overlap=0, p=1.0 | INV-2 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Default background | `FindPathwayEnrichment(query, db)` with no background | background = union(pathways ∪ query); p computed against that N | ASSUMPTION (background defaulting) |
| S2 | Duplicate genes in query | query with repeated gene ids | de-duplicated; QuerySize counts distinct | Set semantics |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | p-value bounded | random small inputs | 0 ≤ p ≤ 1 (INV-1) | Property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `FindPathwayEnrichment` / `HypergeometricUpperTail`. Implementation already
  present in `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`
  (`FindPathwayEnrichment`, `HypergeometricUpperTail`, `LogChoose`, `LogGamma`). No prior tests found.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14, S1–S2, C1 | ❌ Missing | No prior tests existed for this unit |

<!-- Status values: ✅ Covered, ⚠ Weak, ❌ Missing, 🔁 Duplicate -->

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_FindPathwayEnrichment_Tests.cs` — all cases above.
- **Remove:** none (no pre-existing tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| MetagenomicsAnalyzer_FindPathwayEnrichment_Tests.cs | Canonical unit tests | 17 |

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
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented (merged with M8 Assert.Multiple) | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | S1 | ❌ Missing | Implemented | ✅ Done |
| 16 | S2 | ❌ Missing | Implemented | ✅ Done |
| 17 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | PNNL worked example reproduced |
| M2 | ✅ | Symmetric orientation |
| M3 | ✅ | Exact 1/252 |
| M4 | ✅ | Exact 5/6 |
| M5 | ✅ | x=0 → 1 |
| M6 | ✅ | Exact 251/252 |
| M7 | ✅ | Degenerate → 1 (all variants) |
| M8 | ✅ | End-to-end p = 1/12 |
| M9 | ✅ | Ascending sort verified |
| M10 | ✅ | Fields verified (merged into M8) |
| M11 | ✅ | Null query throws |
| M12 | ✅ | Null DB throws |
| M13 | ✅ | Empty DB → empty |
| M14 | ✅ | No overlap → p=1 |
| S1 | ✅ | Default background path |
| S2 | ✅ | Duplicate de-dup |
| C1 | ✅ | Bounded property test |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Default background = union(pathway members ∪ query) when none supplied | S1 |

---

## 7. Open Questions / Decisions

1. Pathway-to-gene mappings are NOT hard-coded (KEGG/MetaCyc are large curated DBs). The unit is a
   generic caller-supplied-definitions ORA calculator; the math is what is tested. Decision recorded.
