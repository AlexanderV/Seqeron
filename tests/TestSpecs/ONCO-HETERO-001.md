# Test Specification: ONCO-HETERO-001

**Test Unit ID:** ONCO-HETERO-001
**Area:** Oncology
**Algorithm:** Tumor Heterogeneity Analysis (MATH score, Shannon clonal diversity, subclone count, subclonal fraction)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Mroz & Rocco (2013), Oral Oncology 49(3):211–215 (MATH) | 1 | https://pubmed.ncbi.nlm.nih.gov/23079694/ | 2026-06-15 |
| 2 | Mroz et al. (2015), PLOS Medicine 12(2):e1001786 (MAD 1.4826) | 1 | https://doi.org/10.1371/journal.pmed.1001786 | 2026-06-15 |
| 3 | maftools `mathScore.R` (reference impl.) | 3 | https://github.com/PoisonAlien/maftools/blob/master/R/mathScore.R | 2026-06-15 |
| 4 | Liu & Zhang (2017), BMC Genomics 18:457 (Shannon ITH) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5468233/ | 2026-06-15 |
| 5 | Shannon (1948), Bell Syst. Tech. J. 27:379–423 | 4 | https://en.wikipedia.org/wiki/Diversity_index#Shannon_index | 2026-06-15 |
| 6 | Landau et al. (2013), Cell 152(4):714–726 (CCF 0.95 threshold) | 1 | https://doi.org/10.1016/j.cell.2013.01.019 | 2026-06-15 |

### 1.2 Key Evidence Points

1. MATH = 100 × MAD/median over mutant-allele fractions — Mroz & Rocco (2013).
2. MAD = 1.4826 × median(|f − median(f)|) (normal consistency) — Mroz et al. (2015); maftools `pat.math = pat.mad * 1.4826 / median(vaf)`.
3. Shannon diversity H = −Σ pᵢ ln(pᵢ), natural log, over clone fractions — Liu & Zhang (2017); Shannon (1948).
4. Richness = number of clones/clusters present — Liu & Zhang (2017).
5. Subclonal mutation ⇔ CCF < 0.95 — Landau et al. (2013).

### 1.3 Documented Corner Cases

- Median MAF = 0 ⇒ MATH undefined (division by zero) — maftools/formula.
- All-identical VAFs or single mutation ⇒ MAD = 0 ⇒ MATH = 0.
- Single clone ⇒ Shannon H = 0; k equal clones ⇒ H = ln k.

### 1.4 Known Failure Modes / Pitfalls

1. Forgetting the 1.4826 scaling underestimates MATH by ~33% — Mroz et al. (2015).
2. Using log2 instead of ln changes H by a constant factor — Shannon (1948); sources use ln.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateITH(ccfDistribution)` | OncologyAnalyzer | Canonical | MATH score = 100·1.4826·MAD/median |
| `InferSubclones(ccfClusters)` | OncologyAnalyzer | Canonical | count of occupied CCF clusters |
| `AnalyzeHeterogeneity(vafs, ccfValues, k)` | OncologyAnalyzer | Canonical | aggregate: MATH + Shannon + subclones + subclonal fraction |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | MATH (ITH_score) ≥ 0 | Yes | MAD ≥ 0 and median > 0 — Mroz & Rocco (2013) |
| INV-2 | MATH = 0 ⇔ MAD = 0 (all values equal the median) | Yes | formula |
| INV-3 | Shannon H ≥ 0; H = 0 ⇔ single occupied clone | Yes | Shannon (1948); Liu & Zhang (2017) |
| INV-4 | Shannon H ≤ ln(richness), equality for equal clones | Yes | Shannon (1948) |
| INV-5 | 1 ≤ subclone count ≤ k; 0 ≤ subclonal fraction ≤ 1 | Yes | Liu & Zhang (2017); Landau (2013) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | MATH odd count | VAFs {0.1,0.2,0.3,0.4,0.5} | 49.42 (=100·1.4826·0.10/0.30) | Mroz & Rocco (2013); maftools |
| M2 | MATH even count | VAFs {0.2,0.4,0.6,0.8} | 59.304 (=100·1.4826·0.20/0.50) | Mroz et al. (2015) |
| M3 | MATH all identical | VAFs {0.3,0.3,0.3} | 0.0 (MAD = 0) | formula / INV-2 |
| M4 | MATH single value | VAFs {0.4} | 0.0 | formula |
| M5 | Shannon two equal clones | 4 CCFs → 2 clusters of size 2 | H = −ln 0.5 = 0.6931471805599453 | Liu & Zhang (2017); Shannon (1948) |
| M6 | Shannon four equal clones | 4 CCFs → 4 clusters size 1 | H = ln 4 = 1.3862943611198906 | Shannon (1948) |
| M7 | Shannon single clone | k = 1 | H = 0.0 | INV-3 |
| M8 | Subclone count | clustering with 3 occupied clusters | 3 | Liu & Zhang (2017) richness |
| M9 | Subclonal fraction | CCFs {0.4,0.5,0.98,1.0}, threshold 0.95 | 0.5 (2 of 4 below 0.95) | Landau et al. (2013) |
| M10 | AnalyzeHeterogeneity aggregate | VAFs+CCFs example, k=2 | MATH, H, subclones, subclonal fraction all match component derivations | sources 1–6 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Zero median throws | VAFs {0,0,0.4} median 0 | ArgumentException | division by zero |
| S2 | Null distribution throws | CalculateITH(null) | ArgumentNullException | guard |
| S3 | Empty distribution throws | CalculateITH(empty) | ArgumentException | guard |
| S4 | Out-of-range VAF throws | VAF 1.5 | ArgumentException | [0,1] domain |
| S5 | Mismatched lengths throw | vafs.Count ≠ ccf.Count | ArgumentException | alignment |
| S6 | InferSubclones empty throws | empty clustering | ArgumentException | guard |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-1 property | MATH ≥ 0 over varied valid inputs | always ≥ 0 | registry invariant |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `CalculateITH`, `InferSubclones`, or `AnalyzeHeterogeneity`. ONCO-HETERO-001 was Not Started; no `OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs` existed.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10 | ❌ Missing | new unit, no prior tests |
| S1–S6 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs` — all cases for the three methods.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_AnalyzeHeterogeneity_Tests.cs | canonical | 17 |

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
| 16 | S6 | ❌ Missing | implemented | ✅ Done |
| 17 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M10 | ✅ Covered | exact evidence-derived values |
| S1–S6 | ✅ Covered | guards/failure modes |
| C1 | ✅ Covered | invariant property test |

Total in-scope cases: 17; ✅ = 17.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Shannon clone fractions = per-cluster mutation proportions | AnalyzeHeterogeneity (Shannon), M5–M7, M10 |
| 2 | Even-count median = mean of two central order statistics (R/maftools) | CalculateITH median, M2 |

---

## 7. Open Questions / Decisions

1. None — all metrics and constants are traceable to retrieved authoritative sources.
