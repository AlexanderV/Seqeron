# Test Specification: ONCO-ASCAT-001

**Test Unit ID:** ONCO-ASCAT-001
**Area:** Oncology
**Algorithm:** Upstream allele-specific derivation — segmentation, joint purity/ploidy fit (ASCAT), mutation multiplicity
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Van Loo et al. (2010), ASCAT, PNAS 107:16910 | 1 | https://doi.org/10.1073/pnas.1009843107 | 2026-06-23 |
| 2 | ascat.runAscat.R (VanLoo-lab/ascat, master) | 3 | https://github.com/VanLoo-lab/ascat | 2026-06-23 |
| 3 | McGranahan et al. (2016), Science 351:1463 | 1 | https://doi.org/10.1126/science.aaf1490 | 2026-06-23 |
| 4 | Zheng et al. (2022), PICTograph, Bioinformatics 38:3677 | 1 | https://doi.org/10.1093/bioinformatics/btac440 | 2026-06-23 |
| 5 | DeCiFering (2021), PMC8542635 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/ | 2026-06-23 |

### 1.2 Key Evidence Points

1. ASCAT raw copy numbers: `nA = (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho`, `nB = (rho-1 + b*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho` — ascat.runAscat.R (source 2).
2. Goodness-of-fit distance = Σ |nMinor − round(nMinor)|² · length · (b==0.5 ? 0.05 : 1); `goodnessOfFit = (1 − d/TheoretMaxdist)·100`, `TheoretMaxdist = Σ 0.25·length·weight` — source 2.
3. Grid search over (ploidy ψ × aberrant fraction ρ), selecting copy numbers "as close as possible to nonnegative whole numbers" — Van Loo 2010 Fig. 1 (source 1).
4. γ = 1 for sequencing data — ASCAT README (source 2).
5. n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)]; CCF = n_mut/M; M (multiplicity) is n_mut rounded for a clonal mutation — McGranahan 2016 (source 3), PICTograph inversion (source 4).

### 1.3 Documented Corner Cases

- Balanced (BAF = 0.5) segments carry little allele-specific information → ×0.05 GoF weight (source 2).
- Multiple sunrise optima (2n vs 4n) — global minimum over the grid is selected (source 1).
- Multiplicity must be clamped to [1, major CN]: an observed variant has ≥ 1 mutated copy (sources 3, 4).

### 1.4 Known Failure Modes / Pitfalls

1. Wrong γ rescales logR and biases copy number — source 2.
2. Symmetric heterozygous BAF clusters (b and 1−b) cancel if averaged naively → mirror about 0.5 before averaging (standard allele-specific summary).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `SegmentAlleleSpecific` | OncologyAnalyzer | **Canonical** | PCF/CBS-style mean-shift segmentation of per-locus logR/BAF |
| `FitPurityPloidy` | OncologyAnalyzer | **Canonical** | ASCAT grid fit; recovers ρ, ψ, integer (nA,nB) |
| `DeriveMultiplicity` | OncologyAnalyzer | **Canonical** | McGranahan n_mut rounding/clamp |
| `EstimateCcf` | OncologyAnalyzer | **Delegate** | already tested in ONCO-CCF-001; here only end-to-end with derived inputs |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | At the planted (ρ₀, ψ₀), the GoF distance is ≈ 0 (integer copy numbers exactly recovered) | Yes | source 1/2 |
| INV-2 | Recovered ρ, ψ equal planted ρ₀, ψ₀ within the grid step | Yes | source 1 |
| INV-3 | Derived multiplicity ∈ [1, majorCopyNumber] | Yes | sources 3, 4 |
| INV-4 | End-to-end CCF of a planted clonal mutation = 1.0 (within tolerance) | Yes | sources 3, 5 |
| INV-5 | Segment count = number of distinct mean-shift runs; breakpoints at planted change positions | Yes | source 1 (segmentation) |
| INV-6 | GoF percentage ≤ 100; distance at true (ρ,ψ) ≤ distance at a wrong (ρ,ψ) | Yes | source 2 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Segmentation recovers breakpoints | Per-locus logR with two clear levels and a chromosome change | 3 segments at the planted boundaries | source 1 |
| M2 | Grid fit recovers ρ₀, ψ₀ | Synthesise (r,b) from ρ₀=0.80, ψ₀=2.2, segs {1+1,2+0,1+1,2+1,1+1} | Purity=0.80, Ploidy=2.2 | source 1/2 |
| M3 | Grid fit recovers integer (nA,nB) | Same input as M2 | Segments: (1,1),(2,0),(1,1),(2,1),(1,1) | source 2 |
| M4 | GoF ≈ 100% at true params | Same input as M2 | GoodnessOfFit ≈ 100 (distance ≈ 0) | source 2 |
| M5 | Multiplicity m=1 recovered | VAF=0.40 from m=1, CN=2(1+1), ρ=0.80 | DeriveMultiplicity = 1 | source 3/4 |
| M6 | Multiplicity m=2 recovered | VAF=4/7 from m=2, CN=3, major=2, ρ=0.80 | DeriveMultiplicity = 2 | source 3/4 |
| M7 | Multiplicity clamped to major CN | High VAF that rounds above major CN | result = majorCopyNumber | source 3/4 |
| M8 | Multiplicity clamped to ≥ 1 | Tiny VAF that rounds to 0 | result = 1 | source 3/4 |
| M9 | End-to-end CCF = 1.0 (clonal) | Fit → derive CN, multiplicity → EstimateCcf on VAF=0.40 | CCF = 1.0 | source 3/5 |
| M10 | DeriveMultiplicity invalid args throw | vaf>1, purity≤0, CN<1, major∉[1,CN] | ArgumentOutOfRangeException | contract |
| M11 | SegmentAlleleSpecific invalid args throw | null loci, threshold≤0, minLoci<1 | ArgumentNullException / ArgumentOutOfRangeException | contract |
| M12 | FitPurityPloidy invalid args throw | null/empty segments, bad grid bounds | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException | contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | GoF discriminates | Distance at true (ρ,ψ) < distance at deliberately wrong (ρ,ψ) | true < wrong | INV-6 |
| S2 | Balanced-only genome | All loci b=0.5 → segments down-weighted ×0.05 | fit completes; balanced segments folded BAF=0.5 | corner case |
| S3 | Triploid planted (ψ₀=3) | Synthesise from ψ₀=3.0 | recovers ψ≈3.0 | aneuploidy |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Single-locus segment | One locus per chromosome | one segment per chromosome, LocusCount=1 | robustness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New unit; no prior tests for `SegmentAlleleSpecific`, `FitPurityPloidy`, `DeriveMultiplicity`. `EstimateCcf` is covered by ONCO-CCF-001; here used only for the end-to-end M9 check.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M12, S1–S3, C1 | ❌ Missing | brand-new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AscatDerivation_Tests.cs`
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_AscatDerivation_Tests.cs | Canonical | 16+ |

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
| 13 | S1 | ❌ Missing | Implemented | ✅ Done |
| 14 | S2 | ❌ Missing | Implemented | ✅ Done |
| 15 | S3 | ❌ Missing | Implemented | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M12 | ✅ | Implemented, evidence-based |
| S1–S3 | ✅ | Implemented |
| C1 | ✅ | Implemented |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Germline-het-SNP BAF forward model b = (ρ·nB + (1−ρ))/(ρ·n + 2(1−ρ)) — algebraic inverse of the ASCAT nA/nB equations | planted-truth synthesis only (M2–M4, M9, S1–S3) |
| 2 | logR baseline = average sample ploidy (segment at genome-average CN ⇒ r = 0) | planted-truth synthesis only |

Both assumptions affect only test-input synthesis (not production code) and are exact inverses of the cited ASCAT equations.

---

## 7. Open Questions / Decisions

1. None. The full pknotsRG-style refinements of ASCAT (ASPCF multi-sample segmentation, sub-clonal copy number, refit heuristics) are out of scope; this unit derives the canonical single-sample (ρ, ψ, integer nA/nB) + multiplicity that the downstream Oncology methods consume.
