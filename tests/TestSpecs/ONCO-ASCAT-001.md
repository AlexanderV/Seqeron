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
| 6 | Nilsen et al. (2012), Copynumber/PCF+ASPCF, BMC Genomics 13:591 | 1 | https://doi.org/10.1186/1471-2164-13-591 | 2026-06-23 |
| 7 | Ross et al. (2021), ASCAT ASPCF, Bioinformatics 37:1909 | 1 | https://doi.org/10.1093/bioinformatics/btaa538 | 2026-06-23 |
| 8 | Nik-Zainal et al. (2012), Battenberg, Cell 149:994 | 1 | https://doi.org/10.1016/j.cell.2012.04.023 | 2026-06-23 |

### 1.2 Key Evidence Points

1. ASCAT raw copy numbers: `nA = (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho`, `nB = (rho-1 + b*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho` — ascat.runAscat.R (source 2).
2. Goodness-of-fit distance = Σ |nMinor − round(nMinor)|² · length · (b==0.5 ? 0.05 : 1); `goodnessOfFit = (1 − d/TheoretMaxdist)·100`, `TheoretMaxdist = Σ 0.25·length·weight` — source 2.
3. Grid search over (ploidy ψ × aberrant fraction ρ), selecting copy numbers "as close as possible to nonnegative whole numbers" — Van Loo 2010 Fig. 1 (source 1).
4. γ = 1 for sequencing data — ASCAT README (source 2).
5. n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)]; CCF = n_mut/M; M (multiplicity) is n_mut rounded for a clonal mutation — McGranahan 2016 (source 3), PICTograph inversion (source 4).
6. PCF penalised least squares: `L(S|y,γ) = Σ_I Σ_j (y_j − ȳ_I)² + γ|S|`; DP recurrence `e_k = min_j (d_jk + e_{j−1} + γ)`, `e_0=0`; default γ=40 — Nilsen 2012 (source 6).
7. ASPCF joint cost `L(S|y₁,y₂,γ) = L(S|y₁,γ) + L(S|y₂,γ)` (common breakpoints, per-track means); BAF mirrored to a single allelic-imbalance track — Nilsen 2012 / Ross 2021 (sources 6, 7).
8. Sub-clonal segment = one (clonal) or two (subclonal) integer states, fractions summing to 1; `n_obs = f·n₁ + (1−f)·n₂` over the bracketing integers — Battenberg (source 8).

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
| `SegmentAlleleSpecificAspcf` | OncologyAnalyzer | **Canonical** | ASPCF penalised-least-squares (PCF DP) joint logR/BAF segmentation [6][7] |
| `FitSubclonalCopyNumber` | OncologyAnalyzer | **Canonical** | Battenberg two-state sub-clonal decomposition [8] |
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
| INV-7 | ASPCF penalised cost is the global minimum: ≤ greedy cost on the same track | Yes | source 6 |
| INV-8 | ASPCF: γ→large ⇒ 1 segment; small γ recovers each level; no segment crosses a chromosome | Yes | source 6 |
| INV-9 | Sub-clonal state fractions sum to 1; integer (nA,nB) ⇒ single clonal state (f=1) | Yes | source 8 |

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
| M-ASPCF-1 | ASPCF recovers single breakpoint | two-level logR (0 then 1, 10+10 loci), γ=0.5 | 2 segments, breakpoint at index 10, means 0.0/1.0 | source 6 |
| M-ASPCF-2 | ASPCF ≤ greedy cost | noisy two-level track | penalised cost(ASPCF) ≤ cost(greedy) | source 6 |
| M-ASPCF-3 | mirrored-BAF splits same-logR | balanced (BAF 0.5) then LOH (BAF 0.0), same logR | 2 segments (mirrored BAF 0.5 vs 1.0) | source 7 |
| M-SUB-1 | sub-clonal mixture recovered | 0.4·(2,0)+0.6·(1,1) observed at ρ=1,ψ=2 | two states (1,1)/(2,0), f≈0.4, fracs sum to 1 | source 8 |
| M-SUB-2 | clonal integer collapses | (nA,nB)=(2,1) integer at ρ=1,ψ=2 | single state (2,1), f=1, no secondary | source 8 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | GoF discriminates | Distance at true (ρ,ψ) < distance at deliberately wrong (ρ,ψ) | true < wrong | INV-6 |
| S2 | Balanced-only genome | All loci b=0.5 → segments down-weighted ×0.05 | fit completes; balanced segments folded BAF=0.5 | corner case |
| S3 | Triploid planted (ψ₀=3) | Synthesise from ψ₀=3.0 | recovers ψ≈3.0 | aneuploidy |
| S-ASPCF-1 | penalty controls segment count | two-level track, γ=1000 vs γ=0.5 | 1 segment vs 2 segments | source 6 |
| S-ASPCF-2 | chromosome boundary not crossed | flat value over chr1+chr2, γ=100 | 2 segments (one per chromosome) | source 6 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Single-locus segment | One locus per chromosome | one segment per chromosome, LocusCount=1 | robustness |
| C-ASPCF-1 | ASPCF invalid args throw | null loci, penalty≤0 | ArgumentNullException / ArgumentOutOfRangeException | contract |
| C-SUB-1 | sub-clonal invalid args throw | null segments, ρ≤0, ψ≤0, γ≤0 | ArgumentNullException / ArgumentOutOfRangeException | contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- New unit; no prior tests for `SegmentAlleleSpecific`, `FitPurityPloidy`, `DeriveMultiplicity`. `EstimateCcf` is covered by ONCO-CCF-001; here used only for the end-to-end M9 check.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M12, S1–S3, C1 | ✅ Covered | prior session (greedy segmentation, fit, multiplicity) |
| M-ASPCF-1..3, S-ASPCF-1..2, C-ASPCF-1 | ❌ Missing | new ASPCF half |
| M-SUB-1..2, C-SUB-1 | ❌ Missing | new sub-clonal half |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AscatDerivation_Tests.cs`
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_AscatDerivation_Tests.cs | Canonical | 23 |

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
| 17 | M-ASPCF-1 | ❌ Missing | Implemented | ✅ Done |
| 18 | M-ASPCF-2 | ❌ Missing | Implemented | ✅ Done |
| 19 | M-ASPCF-3 | ❌ Missing | Implemented | ✅ Done |
| 20 | S-ASPCF-1 | ❌ Missing | Implemented | ✅ Done |
| 21 | S-ASPCF-2 | ❌ Missing | Implemented | ✅ Done |
| 22 | C-ASPCF-1 | ❌ Missing | Implemented | ✅ Done |
| 23 | M-SUB-1 | ❌ Missing | Implemented | ✅ Done |
| 24 | M-SUB-2 | ❌ Missing | Implemented | ✅ Done |
| 25 | C-SUB-1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 25
**✅ Done:** 25 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M12 | ✅ | Implemented, evidence-based |
| S1–S3 | ✅ | Implemented |
| C1 | ✅ | Implemented |
| M-ASPCF-1..3 | ✅ | Implemented, evidence-based (Nilsen 2012 / Ross 2021) |
| S-ASPCF-1..2 | ✅ | Implemented |
| C-ASPCF-1 | ✅ | Implemented |
| M-SUB-1..2 | ✅ | Implemented, evidence-based (Battenberg) |
| C-SUB-1 | ✅ | Implemented |

---

## 6. Assumption Register

**Total assumptions:** 4

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Germline-het-SNP BAF forward model b = (ρ·nB + (1−ρ))/(ρ·n + 2(1−ρ)) — algebraic inverse of the ASCAT nA/nB equations | planted-truth synthesis only (M2–M4, M9, S1–S3, M-SUB-1/2) |
| 2 | logR baseline = average sample ploidy (segment at genome-average CN ⇒ r = 0) | planted-truth synthesis only |
| 3 | ASPCF γ exposed as a sourced parameter (form `+γ\|S\|` verbatim; numeric default copynumber 40 / ASCAT 70 is probe-scale-specific) | `SegmentAlleleSpecificAspcf` default; M-ASPCF tests use a γ derived from each dataset's ΔSSE |
| 4 | Two-state sub-clonal mixture uses the two bracketing integers with one shared fraction (3+-population mixtures out of scope) | `FitSubclonalCopyNumber`, M-SUB-1/2 |

Assumptions 1–2 affect only test-input synthesis (not production code) and are exact inverses of the cited ASCAT equations. Assumptions 3–4 are sourced modelling choices documented in the Evidence (ASPCF γ form is verbatim; the two-state mixture is the unique f∈[0,1] decomposition of a single fractional value).

---

## 7. Open Questions / Decisions

1. None. ASPCF penalised-least-squares segmentation (`SegmentAlleleSpecificAspcf`) and two-state sub-clonal copy number (`FitSubclonalCopyNumber`) are now implemented. Remaining out-of-scope refinements: multi-sample (asmultipcf) segmentation, 3+-population per-segment mixtures, and a whole-genome-doubling refit search — these are documented limitations, not blockers.
