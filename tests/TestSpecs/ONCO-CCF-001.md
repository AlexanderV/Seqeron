# Test Specification: ONCO-CCF-001

**Test Unit ID:** ONCO-CCF-001
**Area:** Oncology
**Algorithm:** Cancer Cell Fraction (CCF) point estimation and 1D CCF clustering into clones/subclones
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Tarabichi et al. 2021, *Nature Methods* (subclonal reconstruction guide, Box 1) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7867630/ | 2026-06-15 |
| 2 | Zheng et al. 2022, *Bioinformatics* (PICTograph) | 1 | https://academic.oup.com/bioinformatics/article/38/15/3677/6596597 | 2026-06-15 |
| 3 | McGranahan et al. 2016, *Science* 351:1463–1469 | 1 | https://www.science.org/doi/10.1126/science.aaf1490 | 2026-06-15 |
| 4 | CNAqc — CCF computation vignette | 3 | https://caravagnalab.github.io/CNAqc/articles/a4_ccf_computation.html | 2026-06-15 |
| 5 | Lloyd 1982, *IEEE Trans. Inf. Theory* 28(2):129–137 (k-means) | 1 | https://doi.org/10.1109/TIT.1982.1056489 | 2026-06-15 |

### 1.2 Key Evidence Points

1. CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m), with VAF=f, ρ=purity, N_T=tumor copy number, normal CN=2, m=multiplicity — Box 1, Tarabichi 2021; corroborated by Zheng 2022 (VAF=(m·CCF·p)/(c·p+2(1−p))) and McGranahan 2016 (n_mut = VAF·(1/p)·[p·CN_t+2(1−p)], CCF = n_mut/m).
2. Multiplicity m = f·(ρ·N_T + 2(1−ρ))/ρ, rounded to the nearest non-zero integer for clonal CN regions — Box 1, Tarabichi 2021.
3. The cluster with the highest cellular prevalence/CCF is the clonal cluster — Tarabichi 2021 SNV clustering section.
4. Raw CCF can exceed 1 due to sampling noise (CNAqc shows 1.06); registry invariant requires 0 ≤ CCF ≤ 1 — CNAqc + registry.
5. Lloyd k-means: assign each point to nearest centroid (least squared distance), recompute centroids as cluster means, minimize WCSS — Lloyd 1982.

### 1.3 Documented Corner Cases

- CCF > 1 from noise (cap at 1). — CNAqc / Tarabichi 2021.
- Multi-copy loci: multiplicity must be supplied (1 ≤ m ≤ tumor copy number). — Tarabichi 2021.
- Unknown/invalid purity: ρ ∈ (0,1] (denominator). — formula domain.

### 1.4 Known Failure Modes / Pitfalls

1. Treating VAF·2 as CCF regardless of copy number/multiplicity overestimates CCF at amplified/LOH loci. — Tarabichi 2021.
2. Non-deterministic k-means seeding gives unstable clusters; fixed quantile seeding required. — Lloyd 1982 (algorithm is seed-dependent).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `EstimateCcf(vaf, purity, tumorCopyNumber, multiplicity)` | OncologyAnalyzer | Canonical | Point CCF per McGranahan/PMC formula; returns raw + capped. |
| `ClusterCcfValues(ccfValues, clusterCount)` | OncologyAnalyzer | Canonical | Deterministic 1D Lloyd k-means; identifies clonal (max-centroid) cluster. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Reported CCF ∈ [0, 1] (raw capped at 1) | Yes | Registry invariant + McGranahan clonal def |
| INV-2 | CCF is monotonically increasing in VAF (other inputs fixed) | Yes | Formula linear in VAF (Tarabichi 2021) |
| INV-3 | Every clustered value is assigned to exactly one cluster; assignments ∈ [0, k) | Yes | Lloyd 1982 partition |
| INV-4 | Cluster centroid equals the mean of its members | Yes | Lloyd 1982 update step |
| INV-5 | Clonal cluster index = argmax centroid | Yes | Tarabichi 2021 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | EstimateCcf case A | f=0.40, ρ=0.80, N_T=2, m=1 | CCF = 1.0 (raw 1.0) | Formula (Tarabichi 2021) |
| M2 | EstimateCcf case B | f=0.20, ρ=0.80, N_T=2, m=1 | CCF = 0.5 | Formula |
| M3 | EstimateCcf case C multi-copy | f=0.50, ρ=1.0, N_T=4, m=2 | CCF = 1.0 | Formula (multi-copy) |
| M4 | EstimateCcf case D purity 0.5 | f=0.25, ρ=0.5, N_T=2, m=1 | CCF = 1.0 | Formula |
| M5 | EstimateCcf raw < 1 | f=0.471, ρ=1.0, N_T=2, m=1 | CCF = 0.942, raw 0.942 | Formula / CNAqc |
| M6 | EstimateCcf cap at 1 | f=0.60, ρ=0.80, N_T=2, m=1 | reported CCF = 1.0, raw = 1.5 | INV-1 / CNAqc |
| M7 | EstimateCcf invalid purity | ρ = 0 or 1.2 | ArgumentOutOfRangeException | Formula domain |
| M8 | EstimateCcf invalid VAF | f = −0.1 or 1.1 | ArgumentOutOfRangeException | VAF ∈ [0,1] |
| M9 | EstimateCcf invalid copy number | N_T = 0 | ArgumentOutOfRangeException | N_T ≥ 1 |
| M10 | EstimateCcf invalid multiplicity | m = 0 or m > N_T | ArgumentException | Multiplicity def |
| M11 | ClusterCcfValues 2 clones | {1.0,0.98,0.96,0.50,0.48,0.52}, k=2 | centroids {0.50,0.98}; low={3,4,5}, high={0,1,2}; clonal=high | Lloyd 1982 + Tarabichi 2021 |
| M12 | ClusterCcfValues determinism | same input shuffled, k=2 | identical centroids and per-value assignment | Deterministic seeding |
| M13 | ClusterCcfValues clonal index | dataset M11 | ClonalClusterIndex centroid = 0.98 (max) | Tarabichi 2021 |
| M14 | ClusterCcfValues null | null input | ArgumentNullException | Failure mode |
| M15 | ClusterCcfValues k invalid | k=0 or k>n | ArgumentOutOfRangeException | k ∈ [1, n] |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | ClusterCcfValues k=1 | {0.3,0.6,0.9}, k=1 | one cluster, centroid=0.6, clonal index 0 | boundary |
| S2 | EstimateCcf monotonicity (INV-2) | f1<f2 same other inputs | CCF(f1) < CCF(f2) | property |
| S3 | ClusterCcfValues empty | empty list | ArgumentException | failure mode |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | EstimateCcf raw zero VAF | f=0 | CCF = 0 | trivial boundary |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing `EstimateCcf` / `ClusterCcfValues` methods or tests in `OncologyAnalyzer.cs` (grep over the class). The neighbouring `ClassifyClonality` (ONCO-CLONAL-001) uses a Bayesian grid posterior — a distinct algorithm; this unit adds the deterministic point estimate and 1D clustering. No duplication.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M15, S1–S3, C1 | ❌ Missing | New unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimateCcf_Tests.cs` — all cases for both methods (single unit).
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_EstimateCcf_Tests.cs` | Canonical (this unit) | 19 |

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
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | M15 | ❌ Missing | Implemented | ✅ Done |
| 16 | S1 | ❌ Missing | Implemented | ✅ Done |
| 17 | S2 | ❌ Missing | Implemented | ✅ Done |
| 18 | S3 | ❌ Missing | Implemented | ✅ Done |
| 19 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 19
**✅ Done:** 19 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact CCF 1.0 |
| M2 | ✅ Covered | exact CCF 0.5 |
| M3 | ✅ Covered | multi-copy CCF 1.0 |
| M4 | ✅ Covered | purity 0.5 CCF 1.0 |
| M5 | ✅ Covered | CCF 0.942 |
| M6 | ✅ Covered | cap + raw |
| M7 | ✅ Covered | purity guard |
| M8 | ✅ Covered | VAF guard |
| M9 | ✅ Covered | copy-number guard |
| M10 | ✅ Covered | multiplicity guard |
| M11 | ✅ Covered | centroids + assignments |
| M12 | ✅ Covered | determinism |
| M13 | ✅ Covered | clonal index |
| M14 | ✅ Covered | null guard |
| M15 | ✅ Covered | k guard |
| S1 | ✅ Covered | k=1 |
| S2 | ✅ Covered | monotonicity |
| S3 | ✅ Covered | empty guard |
| C1 | ✅ Covered | zero VAF |

**In-scope cases:** 19 — **✅:** 19

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Reported CCF capped to [0,1] (raw exposed) per registry invariant + McGranahan clonal def | M6, INV-1 |
| 2 | 1D clustering = deterministic Lloyd k-means with quantile seeding (no RNG); clonal=max centroid | M11–M13, S1, INV-3..5 |

---

## 7. Open Questions / Decisions

1. Multiplicity is an explicit input (caller-supplied integer), matching the multi-region/PICTograph convention; automatic multiplicity inference from VAF is out of scope for this unit (would belong to a copy-number phasing unit).
