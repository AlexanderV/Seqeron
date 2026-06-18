# Test Specification: ONCO-EXPR-001

**Test Unit ID:** ONCO-EXPR-001
**Area:** Oncology
**Algorithm:** Tumor Gene Expression Outlier (z-score) and Signature Score
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | cBioPortal — mRNA z-score normalization spec | 5 | https://docs.cbioportal.org/z-score-normalization-script/ | 2026-06-15 |
| 2 | cBioPortal-core NormalizeExpressionLevels.java (reference impl) | 3 | https://github.com/cBioPortal/cbioportal-core/blob/master/src/main/java/org/mskcc/cbio/portal/scripts/NormalizeExpressionLevels.java | 2026-06-15 |
| 3 | cBioPortal FAQ — default ±2 threshold | 5 | https://docs.cbioportal.org/user-guide/faq/ | 2026-06-15 |
| 4 | Lee et al. (2008) Inferring Pathway Activity — combined z-score | 1 | https://doi.org/10.1371/journal.pcbi.1000217 | 2026-06-15 |
| 5 | GSVA vignette — corroborates combined z-score | 3 | https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html | 2026-06-15 |

### 1.2 Key Evidence Points

1. Per-gene expression z-score: z = (r − μ)/σ, μ and σ over the reference cohort — cBioPortal spec (source 1) and `getZ` in NormalizeExpressionLevels.java (source 2).
2. σ is the **sample** standard deviation with divisor (n−1) — NormalizeExpressionLevels.java `std()`: `std/(double)(v.length-1)` (source 2). This resolves the n vs n−1 ambiguity the prose leaves open.
3. Outlier classification is strict: z > 2 ⇒ overexpressed, z < −2 ⇒ underexpressed; ±2 inclusive is NOT altered — cBioPortal FAQ (source 3).
4. Zero-SD reference is a failure mode: reference impl aborts (`fatalError "… standard deviation of 0.0"`) (source 2).
5. Signature / pathway activity combined z-score: a = (Σᵢ zᵢ)/√k over the k member genes — Lee et al. (2008) (source 4), corroborated by GSVA (source 5).

### 1.3 Documented Corner Cases

- Zero-SD (constant) reference cohort → no defined z-score (reference impl aborts). [source 2]
- Reference cohort of size ≤ 1 → sample SD (n−1) undefined. [source 2]
- z = ±2 boundary → not an outlier (rule is strict). [source 3]
- Single-gene signature (k=1) → a = z₁; empty signature (k=0) → a undefined. [source 4]

### 1.4 Known Failure Modes / Pitfalls

1. Using population SD (n) instead of sample SD (n−1) gives wrong z-scores — reference impl uses (n−1). [source 2]
2. Treating ±2 as altered (≥ instead of >) over-calls outliers. [source 3]
3. Omitting the √k denominator in the signature score (plain mean would divide by k) inflates/deflates the activity. [source 4]

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateExpressionZScore(double, IReadOnlyList<double>)` | OncologyAnalyzer | Canonical | z=(x−μ)/σ, σ sample SD (n−1) |
| `IdentifyOutlierGenes(IReadOnlyDictionary<string,double>, IReadOnlyDictionary<string,IReadOnlyList<double>>, double)` | OncologyAnalyzer | Canonical | per-gene z + strict >thr / <−thr classification |
| `CalculateSignatureScore(IReadOnlyList<double>)` | OncologyAnalyzer | Canonical | combined z-score a = Σz/√k |
| `ExpressionDirection` (enum), `ExpressionOutlier` (record) | OncologyAnalyzer | Internal | result types, tested via IdentifyOutlierGenes |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | z(μ) = 0 when value equals the cohort mean | Yes | z=(μ−μ)/σ=0 [source 1] |
| INV-2 | z is monotone increasing in the value for a fixed cohort | Yes | z linear in x, σ>0 [source 1] |
| INV-3 | z(2μ−x) = −z(x): reflecting the value about the mean negates z | Yes | linearity of z [source 1] |
| INV-4 | Outlier iff z > +threshold (Over) or z < −threshold (Under); strict inequality | Yes | cBioPortal FAQ ">2 or <-2" [source 3] |
| INV-5 | Signature score of k equal z-scores all = c equals c·√k | Yes | a = (k·c)/√k = c√k [source 4] |
| INV-6 | Signature score with k=1 equals the single z-score | Yes | a = z/√1 = z [source 4] |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | z over-expressed | cohort {2,2,4,6,6} (μ=4,σ=2), x=10 | z = 3.0 | source 1,2 |
| M2 | z at mean | same cohort, x=4 | z = 0.0 | source 1 (INV-1) |
| M3 | z under-expressed | same cohort, x=−1 | z = −2.5 | source 1,2 |
| M4 | sample SD (n−1) used | cohort {2,2,4,6,6}: σ=√(16/4)=2, not √(16/5)=1.788 | x=6 ⇒ z=1.0 (n−1), would be 1.118 if n used | source 2 |
| M5 | outlier classify — over | IdentifyOutlierGenes, gene z=3.0, thr=2 | outlier, Over | source 3 |
| M6 | outlier classify — under | gene z=−2.5, thr=2 | outlier, Under | source 3 |
| M7 | boundary not outlier | x=8 ⇒ z=2.0 exactly, thr=2 | NOT an outlier | source 3 (strict >) |
| M8 | non-outlier excluded | x=4 ⇒ z=0 | not in outlier list | source 3 |
| M9 | signature combined z-score | z={3,1,−1,1}, k=4 | a = 4/√4 = 2.0 | source 4 |
| M10 | single-gene signature | z={2.5}, k=1 | a = 2.5 | source 4 (INV-6) |
| M11 | zero-SD reference throws | cohort {5,5,5} | ArgumentException | source 2 (fatalError) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | null reference cohort | CalculateExpressionZScore(x, null) | ArgumentNullException | guard |
| S2 | reference size 1 | cohort {5} | ArgumentException | (n−1) undefined |
| S3 | empty signature | CalculateSignatureScore({}) | ArgumentException | √k=0, k=0 invalid [source 4] |
| S4 | null signature | CalculateSignatureScore(null) | ArgumentNullException | guard |
| S5 | IdentifyOutlierGenes null args | null sample / null cohorts | ArgumentNullException | guard |
| S6 | gene missing reference | sample gene absent from cohorts dict | ArgumentException | cannot z-score without reference |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-3 reflection | z(2μ−x) = −z(x) for cohort {2,2,4,6,6}, x=10 → 2μ−x=−2 → z=−3 | property | symmetry |
| C2 | INV-5 equal-z signature | z={c,c,c,c}=√k·c, c=1.5,k=4 → a=3.0 | property | variance stabilization |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing expression / outlier / signature-score code in `OncologyAnalyzer.cs` (grep for expression/outlier/signature/z-score found only unrelated VAF/CI and SBS-signature code). No prior tests for ONCO-EXPR-001. New unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M11 | ❌ Missing | new unit, no tests exist |
| S1–S6 | ❌ Missing | new unit |
| C1–C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_IdentifyOutlierGenes_Tests.cs` — all cases for the three canonical methods.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_IdentifyOutlierGenes_Tests.cs | canonical | 19 |

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
| 11 | M11 | ❌ Missing | implemented | ✅ Done |
| 12 | S1 | ❌ Missing | implemented | ✅ Done |
| 13 | S2 | ❌ Missing | implemented | ✅ Done |
| 14 | S3 | ❌ Missing | implemented | ✅ Done |
| 15 | S4 | ❌ Missing | implemented | ✅ Done |
| 16 | S5 | ❌ Missing | implemented | ✅ Done |
| 17 | S6 | ❌ Missing | implemented | ✅ Done |
| 18 | C1 | ❌ Missing | implemented | ✅ Done |
| 19 | C2 | ❌ Missing | implemented | ✅ Done |

**Total items:** 19
**✅ Done:** 19 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M11 | ✅ | implemented with exact evidence-derived values |
| S1–S6 | ✅ | argument validation covered |
| C1–C2 | ✅ | property invariants covered |

In-scope cases: 19; ✅ = 19.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Caller supplies reference cohorts and signature gene sets (scope decision, not numeric) | API shape |
| 2 | Inputs are on a normalization scale where a z-score is meaningful (scale-agnostic; no output change) | Contract §3 |

---

## 7. Open Questions / Decisions

1. SD divisor (n vs n−1) — RESOLVED: reference implementation uses (n−1) sample SD.
2. Threshold strictness — RESOLVED: strict `>`/`<` per cBioPortal FAQ wording.
