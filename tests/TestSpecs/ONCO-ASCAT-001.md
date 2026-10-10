# Test Specification: ONCO-ASCAT-001

**Test Unit ID:** ONCO-ASCAT-001
**Area:** Oncology
**Algorithm:** Upstream allele-specific derivation — segmentation, joint purity/ploidy fit (ASCAT), mutation multiplicity
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-28

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
2. Goodness-of-fit distance = Σ |nMinor − max(round(nMinor),0)|² · length · (b==0.5 ? 0.05 : 1) over **autosomal** segments, `length` = number of heterozygous probes (`make_segments`), nMinor chosen genome-wide (`sum(nA) < sum(nB)`); `goodnessOfFit = (1 − d/TheoretMaxdist)·100`, `TheoretMaxdist = Σ 0.25·length·weight` — source 2.
3. Grid ψ ∈ seq(min_ploidy−0.5, max_ploidy+0.5, 0.05) × ρ ∈ seq(0.1, 1.05, 0.01); candidates = strict minima of a 7×7 window; four-pass filter cascade (ploidy range, ρ ≥ 0.2, GoF > 80, percentzero > 0.02 / perczeroAbb > 0.1 / percOddEven > 0.05, non-aberrant flag, ρ>1 columns masked in pass 3); smallest distance wins; ρ>1 reported as 1; no candidate ⇒ rho = NA — `runASCAT` (source 2).
4. γ = 1 for sequencing data — ASCAT README (source 2).
5. n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)]; CCF = n_mut/M; M (multiplicity) is n_mut rounded for a clonal mutation — McGranahan 2016 (source 3), PICTograph inversion (source 4).
6. PCF penalised least squares: `L(S|y,γ) = Σ_I Σ_j (y_j − ȳ_I)² + γ|S|`; DP recurrence `e_k = min_j (d_jk + e_{j−1} + γ)`, `e_0=0` — Nilsen 2012 (source 6). ASCAT (`ascat.aspcf`/`fastAspcf`/`aspcfpart`) standardises each track by its MAD sd (`getMad`), uses kmin = 6, MAD winsorisation (τ 2.5, k 25), 1000-locus windows, BAF level 0.5+μ with shrinkage, penalty 70 (source 2).
7. ASPCF joint cost `Σ SSE₁/sd₁² + SSE₂/sd₂² + γ·#breakpoints` (common breakpoints, per-track means); BAF mirrored to a single allelic-imbalance track — Nilsen 2012 / Ross 2021 / ascat.aspcf.R (sources 6, 7, 2).
8. Sub-clonal segment = one (clonal) or two (subclonal) integer states on the nearest **edge** of the copy-number square (states differ in one allele, `orderEdges` option 1), fraction τ solving the BAF mixture; clonal iff |l − nearest-corner BAF| < maxdist 0.01 or the per-SNP t-test is not significant (constant BAF ⇒ pval 0) — Battenberg `determine_copynumber` (source 8).

### 1.3 Documented Corner Cases

- Balanced (BAF = 0.5) segments carry little allele-specific information → ×0.05 GoF weight (source 2).
- Multiple sunrise optima (2n vs 4n) — ASCAT takes strict local minima that pass the filter cascade and keeps the smallest distance; no acceptable minimum ⇒ rho = NA (`TryFitPurityPloidy` false / `FitPurityPloidy` throws) (source 2).
- Sex chromosomes (X, Y) are excluded from the fit (source 2).
- Multiplicity must be clamped to [1, major CN]: an observed variant has ≥ 1 mutated copy (sources 3, 4).

### 1.4 Known Failure Modes / Pitfalls

1. Wrong γ rescales logR and biases copy number — source 2.
2. Symmetric heterozygous BAF clusters (b and 1−b) cancel if averaged naively → mirror about 0.5 before averaging (standard allele-specific summary).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `SegmentAlleleSpecific` | OncologyAnalyzer | **Delegate** | ≡ `SegmentAlleleSpecificAspcf(loci, 70)` since B24 F35 (former unsourced greedy heuristic removed; thresholds ignored) |
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
| INV-5 | Segmentation = ascat.aspcf output (breakpoints at the R-verified positions); `SegmentAlleleSpecific` ≡ ASPCF at γ = 70 (F35) | Yes | source 2 (R run) |
| INV-6 | GoF percentage ≤ 100; distance at true (ρ,ψ) ≤ distance at a wrong (ρ,ψ) | Yes | source 2 |
| INV-7 | ASPCF penalised cost is the global minimum over segmentations of the track (pinned by the R-locked outputs) | Yes | source 6 |
| INV-8 | ASPCF: γ→large ⇒ 1 segment; small γ recovers each level; no segment crosses a chromosome | Yes | source 6 |
| INV-9 | Sub-clonal state fractions sum to 1; integer (nA,nB) ⇒ single clonal state (f=1) | Yes | source 8 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | `SegmentAlleleSpecific` = ascat.aspcf (penalty 70) (F35) | M-ASPCF-1 noisy step track, logRChangeThreshold 0.5 | 2 segments (40/40), logR −0.023564999999999999/0.60321000000000002, BAF 0.5/0.75129407874999998; equal to `SegmentAlleleSpecificAspcf(…, 70)` | source 2 (R run) |
| M1b | Legacy thresholds ignored (F35) | same track, thresholds (0.05, 0.01, 1) vs (10, 5, 50) | identical segmentations | F35 delegation |
| M2 | Grid fit recovers ρ₀, ψ₀ | Synthesise (r,b) from ρ₀=0.80, ψ₀=2.2, segs {1+1,2+0,1+1,2+1,1+1} | Purity=0.80, Ploidy=2.2 | source 1/2 |
| M3 | Grid fit recovers integer (nA,nB) | Same input as M2 | Segments: (1,1),(2,0),(1,1),(2,1),(1,1) | source 2 |
| M4 | GoF ≈ 100% at true params | Same input as M2 | GoodnessOfFit ≈ 100 (distance ≈ 0) | source 2 |
| M5 | Multiplicity m=1 recovered | VAF=0.40 from m=1, CN=2(1+1), ρ=0.80 | DeriveMultiplicity = 1 | source 3/4 |
| M6 | Multiplicity m=2 recovered | VAF=4/7 from m=2, CN=3, major=2, ρ=0.80 | DeriveMultiplicity = 2 | source 3/4 |
| M7 | Multiplicity clamped to major CN | High VAF that rounds above major CN | result = majorCopyNumber | source 3/4 |
| M8 | Multiplicity clamped to ≥ 1 | Tiny VAF that rounds to 0 | result = 1 | source 3/4 |
| M9 | End-to-end CCF = 1.0 (clonal) | Fit → derive CN, multiplicity → EstimateCcf on VAF=0.40 | CCF = 1.0 | source 3/5 |
| M11 | Multiplicity exact .5 ties → half-to-even (FIN-B24 F26) | (VAF,ρ,N_T,major) = (0.625,1,4,4); (0.375,0.5,2,2); (0.875,1,4,4); (0.5625,1,8,8); (0.125,1,4,4) | 2; 2; 4; 4; 1 (R `expected_mutant_copies` output) | facets-suite `ccf-annotate-maf.R` |
| M10 | DeriveMultiplicity invalid args throw | vaf>1, purity≤0, CN<1, major∉[1,CN] | ArgumentOutOfRangeException | contract |
| M11 | SegmentAlleleSpecific invalid args throw | null loci, threshold≤0, minLoci<1; NaN/+∞ logR, BAF −0.3 / 1.5 (F35) | ArgumentNullException / ArgumentOutOfRangeException / ArgumentException | contract; ASPCF input contract |
| M12 | FitPurityPloidy invalid args throw | null/empty segments, bad grid bounds | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException | contract |
| M-ASPCF-1 | ASPCF = ascat.aspcf on a noisy step | 80 loci, logR 0→0.6, BAF balanced→0.5±0.25, penalty 70 | 2 segments (40/40), logR −0.023565/0.60321, BAF 0.5 (shrunk)/0.751294079 | source 2 (R run) |
| M-ASPCF-2 | noise-free track | 10×0.0 + 10×1.0, penalty 0.5 | 1 segment (MAD 0 ⇒ window skipped), logR 0.5 | source 2 (R run) |
| M-ASPCF-3 | mirrored-BAF splits same-logR | balanced then LOH (0.5±0.47), identical logR, penalty 70 | 2 segments, BAF 0.5 / 0.966975 | source 2 (R run) |
| M-ASCAT-1..2 | FitPurityPloidy = runASCAT | two noisy genomes (one with chrX) | ρ/ψ/GoF/segments of runASCAT (1/2.7/99.781420571107006/2:1×3; 0.85/2.2/99.999772627448223/2:0 2:1 1:1) | source 2 (R run) |
| M-ASCAT-3 | no ASCAT optimum | two genomes where runASCAT returns NA | Try = false; Fit throws InvalidOperationException | source 2 (R run) |
| M-ASCAT-4 | seg_raw rounding | balanced odd total 3; negative allele | 2:1; 2:0 | source 2 |
| M-SUB-1 | edge mixture recovered | 0.3·(2,1)+0.7·(1,1) at ρ=0.8, ψ=2.5 | (1,1)@0.70000000000000051, (2,1)@0.29999999999999949 | source 8 (R run) |
| M-SUB-2 | both alleles fractional | ρ=1, ψ=2, logR 0, BAF 0.7 | (2,0)@0.14285714285714241, (2,1) | source 8 (R run) |
| M-SUB-3 | clonal corners | forward (2,1)/(2,0)/(3,1) at ρ=0.8, ψ=2.5 | single state, f=1 | source 8 (R run) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | GoF discriminates | `EvaluatePurityPloidy` (ASCAT rho_manual/psi_manual) at true vs wrong (ρ,ψ) | 100 % > wrong | INV-6 |
| S2 | Balanced-only genome | All loci b=0.5 → segments down-weighted ×0.05 | fit completes; balanced segments folded BAF=0.5 | corner case |
| S3 | Triploid planted (ψ₀=3) | Synthesise from ψ₀=3.0 | recovers ψ≈3.0 | aneuploidy |
| S-ASPCF-1 | penalty controls segment count | noisy step, γ=1e6 vs default 70 | 1 segment (logR 0.2898225, BAF 0.637905789375) vs 2 | source 2 (R run) |
| S-ASPCF-3 | < 6 loci on a chromosome | 5 loci | 1 segment, logR 0.13, BAF 0.774 (mean mirrored, no shrink) | source 2 (R run) |
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
| M1–M12, S1–S3, C1 | ✅ Covered | prior session (greedy segmentation, fit, multiplicity); M1/M11 re-locked to ASPCF by F35 (§9) |
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
| 3 | (superseded 2026-09) ASPCF penalty is on ASCAT's MAD-standardised cost; default 70 = `ascat.aspcf` | `SegmentAlleleSpecificAspcf` |
| 4 | (superseded 2026-09) Sub-clonal states follow Battenberg's nearest edge; a segment summary has constant SNP BAF ⇒ Battenberg pval = 0, so only the maxdist rule decides clonality | `FitSubclonalCopyNumber` |

Assumptions 1–2 affect only test-input synthesis (not production code) and are exact inverses of the cited ASCAT equations. Assumptions 3–4 are sourced modelling choices documented in the Evidence (ASPCF γ form is verbatim; the two-state mixture is the unique f∈[0,1] decomposition of a single fractional value).

---

## 7. Open Questions / Decisions

1. None. ASPCF penalised-least-squares segmentation (`SegmentAlleleSpecificAspcf`) and two-state sub-clonal copy number (`FitSubclonalCopyNumber`) are now implemented. Remaining out-of-scope refinements: multi-sample (asmultipcf) segmentation, 3+-population per-segment mixtures, and a whole-genome-doubling refit search — these are documented limitations, not blockers.

---

## 8. 2026-09 review (B24, F12–F14)

- **F12** `FitPurityPloidy` is now a line-by-line port of `runASCAT` (distance matrix, probe-count weights, autosomes only,
  7×7 strict local minima, 4-pass filter cascade, ρ>1 ⇒ 1, `seg_raw` rounding with negative correction and the balanced
  odd-total rule, R half-to-even rounding). Additive API: `TryFitPurityPloidy`, `EvaluatePurityPloidy`
  (rho_manual/psi_manual), `PurityPloidyFit.Psi`, `PurityPloidyFit.IsNonAberrant`; `Ploidy` is ASCAT's output ploidy.
  Defaults = ASCAT (purity 0.1–1.05, ploidy 1.5–5.5). R cross-check: 150/150 random genomes identical (pre-fix code
  disagreed on 123/148).
- **F13** `SegmentAlleleSpecificAspcf` is a port of `ascat.aspcf`/`fastAspcf`/`aspcfpart` (+ `madWins`, `getMad`,
  R `runmed`/`smoothEnds`). R cross-check: 60/60 genomes, 681 segments identical to ≤ 5e-16 (pre-fix: 52/60 differ).
- **F14** `FitSubclonalCopyNumber` is a port of Battenberg `determine_copynumber` + `orderEdges`. R cross-check:
  402/402 segments identical.
- Tests: unit (`OncologyAnalyzer_AscatDerivation_Tests.cs`) M-ASCAT-1..4, M-ASPCF-1..3, S-ASPCF-1..4, M-SUB-1..3;
  fuzz (`OncologyAscatFuzzTests.cs`) single-segment genomes ⇒ ASCAT NA, random genomes ⇒ well-formed or NA;
  properties (`OncologyProperties.cs`) ranges/determinism over found optima; metamorphic ASPCF logR-shift on a noisy track.

## 9. FIN-B24 F35 — `SegmentAlleleSpecific` delegates to ASPCF (2026-10-09)

- The unsourced greedy mean-shift heuristic behind `SegmentAlleleSpecific` was removed; the method now returns
  `SegmentAlleleSpecificAspcf(loci, AspcfDefaultPenalty = 70)` (ascat.aspcf port, R-verified F13). Legacy parameters
  `logRChangeThreshold` / `bafChangeThreshold` / `minLociPerSegment` are ignored (still range-checked). Loci with
  NaN/±∞ logR or BAF ∉ [0, 1] → `ArgumentException` (ASPCF input contract).
- M1 now locks the R ascat.aspcf output of the M-ASPCF-1 track through `SegmentAlleleSpecific`; M1b: thresholds
  ignored; M11 extended with the ASPCF validation.
- Fit tests (M2–M4, M9, S1–S3) and the Algebraic / Combinatorial / Metamorphic allele-swap fixtures place each planted
  segment on its own chromosome (5 loci < kmin 6 ⇒ one ascat.aspcf segment per chromosome, S-ASPCF-3; noise-free
  data has MAD 0 ⇒ no in-chromosome breakpoint, S-ASPCF-2) — the planted fit expectations are unchanged.
- Heavy-tier contracts changed: Fuzz `SegmentAlleleSpecific_LohVsBalancedSameLogR_SplitsOnBaf` (now the R-locked
  80-locus LOH track ⇒ 2 segments, 40/40), `SegmentAlleleSpecific_RandomSignal_PreservesLocusCountNoMalformedSummary`
  (out-of-domain BAF on every 5th seed ⇒ `ArgumentException`; in-domain ⇒ invariants + BAF ∈ [0.5, 1]);
  Metamorphic `Ascat_ConstantLogRShift_PreservesGreedyBreakpoints` → `…PreservesSegmentAlleleSpecificBreakpoints`
  (noisy track, R-verified shifts).

## 10. FIN-B24 F36 — `ascat.aspcf` with germline genotypes (2026-10-10)

Test file: `Unit/Oncology/OncologyAnalyzer_AscatGermlineHomozygous_Tests.cs`. Inputs: genomes G1–G3 from the shared
MINSTD generator (Evidence § F36); expected values = output of the original R `ascat.aspcf` (R 4.3.3).

| ID | Test | Input | Expected | Evidence |
|----|------|-------|----------|----------|
| M-ASPCF-G1 | homozygous deletion inside a homozygous stretch | G1 (1080 probes, 608 het) | 11 logR segments = R (extents, levels ≤ 1e-14, first-het BAF ≤ 1e-15); deletion 281000–330000 level −1.9130208073104829, no BAF; chr3 (no het) one segment, no BAF | source 2 (R run) |
| M-ASPCF-G2 | CN breakpoint inside a homozygous gap | G2 (1002 probes) | 7 segments = R; chr1 split at probe 240; homozygous-only chr4 no BAF; 12-probe chr3 | source 2 (R run) |
| M-ASPCF-G3 | focal gain inside a homozygous stretch | G3 (1100 probes) | 8 segments = R; 481000–540000 level 0.26558251092493651, no BAF | source 2 (R run) |
| S-ASPCF-G1 | per-locus tracks | G1 | segments tile all 1080 loci; BAF NaN exactly at homozygous loci; every locus carries its segment level | ascat.aspcf output contract |
| S-ASPCF-G2 | all heterozygous ⇒ het-only overload | het subsets of G1–G3 | segments bit-identical to `SegmentAlleleSpecificAspcf(loci, 70)` | F13 behaviour preserved |
| S-ASPCF-G3 | homozygous probes carry CN information | G1 with homozygous probes dropped | het-only overload has no level < −1; germline overload finds the deletion at 281000 | source 2 (R run) |
| C-ASPCF-G1 | invalid arguments | null loci / genotypes, length mismatch, het BAF NaN, +∞ logR, penalty 0; homozygous BAF NaN accepted | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException | contract |

## 11. FIN-B24 F37 — runASCAT with homozygous segments (2026-10-10)

Same file and genomes as §10; expected values = original R `ascat.runAscat(gamma = 1)` on the `ascat.aspcf` output.

| ID | Test | Input | Expected | Evidence |
|----|------|-------|----------|----------|
| M-ASCAT-G1 | end-to-end fit, homozygous deletion + no-het chromosome | G1 | ρ 0.7, ψ 2.45, ploidy 2.4722222222222223, GoF 99.937918942624279, seg_raw 1:1 2:1 2:1 0:0 2:1 2:1 1:0 2:2 2:0 3:1 1:1 | source 2 (R run) |
| M-ASCAT-G2 | breakpoint in homozygous gap, chrX | G2 | ρ 0.53, ψ 2.35, ploidy 2.1816367265469063, GoF 99.634551649302338, 1:1 3:0 1:0 2:1 2:2 1:0 1:1 | source 2 (R run) |
| M-ASCAT-G3 | focal gain in homozygous stretch | G3 | ρ 0.92, ψ 3.15, ploidy 3.3818181818181818, GoF 98.658353206884712, 2:2 3:1 2:0 4:0 2:0 3:2 2:1 3:0 | source 2 (R run) |
| S-ASCAT-G1 | segment without het probes | G1 | minor 0 on every no-BAF segment; chr3 2:0 (`bafke = 0`) | source 2 |
| S-ASCAT-G2 | ploidy over all probes | G1–G3 | = LocusCount-weighted mean; het-probe mean = R 2.6134868421052633 / 2.2716666666666665 / 3.539047619047619 (≠ ploidy) | source 2 (R run) |
| S-ASCAT-G3 | manual (ρ, ψ) at the optimum | G3 | `EvaluatePurityPloidyFromAspcf` = fit (GoF, segments, ploidy); `TryFit…` true | rho_manual path |
| S-ASCAT-G4 | all heterozygous ⇒ summary fit | het subsets of G1–G3 | identical ρ, ψ, GoF, ploidy, segments to `FitPurityPloidy(SegmentAlleleSpecificAspcf(…))` | F12 behaviour preserved |
| C-ASCAT-G1 | invalid arguments | null; homozygous-only; X-only; purityStep 0; ρ 1.5 | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException | contract |

## 12. FIN-B24 F38 — runASCAT sex-chromosome model: male haploid X/Y (2026-10-10)

File `OncologyAnalyzer_AscatSexChromosome_Tests.cs`; genomes M1–M3 = G1–G3 autosomes + X (PAR / non-PAR) and Y blocks
(generator of §10 plus per-block positions). Expected values = original R `runASCAT` (gamma = 1) with
`ASCATobj$gender` / `ASCATobj$X_nonPAR` set after `ascat.aspcf(X_nonPAR = NULL)`. Models: M1 `Male` (X_nonPAR NULL),
M2 XY with custom non-PAR [890000, 1090000], M3 `MaleWithXNonPar(GRCh37)` (hg19 constants).

| ID | Test | Input | Expected | Evidence |
|----|------|-------|----------|----------|
| M-ASCAT-X1 | germline-aware fit, XY | M1 `Male` | ρ 0.7, ψ 2.45, ploidy 2.181159420289855, GoF 99.937918942624279; X 1:0 2:0 2:0, Y 1:0 1:0 0:0 | source 2 (R run) |
| M-ASCAT-X2 | germline-aware fit, XY, custom non-PAR | M2 | ρ 0.53, ψ 2.35, ploidy 2.0981228668941978; X 3:0 3:0 1:0 2:0 2:1 (PAR1 segment 17001/34001 > ½ ⇒ haploid; PAR2 12001/45001 ⇒ diploid), Y 1:0 | source 2 (R run) |
| M-ASCAT-X3 | germline-aware fit, XY, hg19 non-PAR | M3 | ρ 0.92, ψ 3.15, ploidy 2.9482758620689653; X 1:1 1:1 1:0 3:0 1:1, Y 2:0 0:0 | source 2 (R run) |
| M-ASCAT-X4..X6 | same genomes, XX | M1–M3 `Female` | ploidy 2.1775362318840581 / 1.8805460750853242 / 2.9482758620689653; ρ, ψ, GoF as XY | source 2 (R run, gender "XX") |
| M-ASCAT-X7..X8 | XY without X_nonPAR (whole X haploid) | M2, M3 `Male` | ploidy 2.1365187713310578 (PAR2 4:0) / 2.9482758620689653 (PAR 2:0) | source 2 (R run) |
| M-ASCAT-X9..X14 | summary path (`FitPurityPloidy`), XY and XX | het-only M1–M3 | M1 2.4579831932773111 (XY = XX), M2 2.3914529914529914 / 2.3384615384615386, M3 ρ 0.91 ψ 3.2 3.4351687388987568 (X 100000–155135000 2:0 vs 1:1) | source 2 (R run) |
| S-ASCAT-X1 | gender-less = `Female` | M1–M3, all six entry points | bit-identical fit and segments | default unchanged |
| S-ASCAT-X2 | sex model does not change the fit | M1–M3 | ρ, ψ, GoF, non-aberrant and autosomal segments identical; X/Y nB = 0 | `autoprobes` |
| S-ASCAT-X3 | manual (ρ, ψ) with a sex model | M2, het-only M3 | = fit; `TryFit…` overloads true | rho_manual path |
| S-ASCAT-X4 | hand-computed equations + 50 % rule | ρ 0.5, ψ 2, r 0 | haploid 3:0, diploid 1:1; overlap 5/10 ⇒ diploid, 6/10 ⇒ haploid; Y always haploid; XX ignores X_nonPAR; ploidy 8/3 | `diploidprobes_fixnonPAR` |
| S-ASCAT-X5 | haploid deletion | Y r = −5 | 0:0 (nAraw −0.875) | negative-value correction |
| S-ASCAT-X6 | presets | — | hg19 (2699521, 154931043), hg38 (2781480, 155701382), CHM13 (2394411, 153925834) | ascat.loadData.R |
| C-ASCAT-X1 | invalid arguments | Start > End; undefined gender / genome; null model on all six overloads | ArgumentOutOfRangeException / ArgumentNullException | contract |

## 13. FIN-B24 F39 — Battenberg per-SNP clonality t-test (2026-10-10)

File `OncologyAnalyzer_AscatDerivation_Tests.cs` (region `FitSubclonalCopyNumberWithSnpTest`) and
`Unit/Core/StatisticsHelper_StudentT_Tests.cs`. Expected values = Battenberg `determine_copynumber` sourced verbatim
(R 4.3.3; maxdist 0.01, siglevel 0.05) and R `pbeta` / `pt` / `t.test` (Evidence § F39).

| ID | Test | Input | Expected | Evidence |
|----|------|-------|----------|----------|
| M-SUBT-1 | Battenberg rows | 18 segments of 4 genomes (G1S28/40–45, G2S15/41/42, G3S31/40–42, G4S27/41/42/44) | R pval (≤ 1e−12 rel), clonal/sub-clonal call, states, fractions | R run |
| M-SUBT-2 | t-test overrides maxdist | G1S41 | maxdist-only sub-clonal, Battenberg clonal 2:1 (p 0.5609) | R run |
| M-SUBT-3 | siglevel / maxdist parameters | G1S28 siglevel 0.1; G1S45 maxdist 0.001 / 0 | (2,0)@0.61666732876527997 + (2,1); p 0.46228196278297101 clonal | R run |
| M-SUBT-4 | SNP BAFs not mirrored | BAFseg < 0.5, SNPs 1 − G1S42 | p 1.7441626943676701e−09, (2,0)@0.15605024544591201 | R run |
| S-SUBT-5 | constant / single SNP ≡ maxdist-only overload | 3 segments | fits identical | `sd == 0 ⇒ pval 0` |
| C-SUBT-1 | invalid arguments | null list / SNP list, SNP BAF 1.2 / NaN, ρ 0, siglevel 1.5, maxdist −0.1 | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException | contract |
| M-STAT-T1 | `RegularizedIncompleteBeta` | 9 (x, a, b) incl. a = b = 10⁴, a = 500 b = ½ | R `pbeta` ≤ 1e−13 rel | R run |
| M-STAT-T2 | `StudentTCdf` | 11 (t, ν) incl. ν 2·10⁵, P 1.7e−24 | R `pt` ≤ 1e−13 rel | R run |
| M-STAT-T3 | `OneSampleTTestPValue` | 5 samples | R `t.test(x, mu)$p.value` ≤ 1e−12 rel; n < 2 / constant ⇒ NaN (R stops) | R run |

## 14. FIN-B24 F40 — Battenberg phased-BAF segmentation (2026-10-10)

File `OncologyAnalyzer_BattenbergPhasedSegmentation_Tests.cs` and `Unit/Core/StatisticsHelper_QuantileType7_Tests.cs`.
Expected values = Battenberg `segment.baf.phased` / `determine_copynumber` sourced verbatim (R 4.3.3), inputs from a
deterministic LCG generator reproduced in the test (Evidence § F40).

| ID | Test | Input | Expected | Evidence |
|----|------|-------|----------|----------|
| M-PHSEG-1 | `runFastPcf` < 1000 SNPs | t1 (300 SNPs, 15-SNP switched blocks) | 3 segments, BAFseg ≤ 1e−12, Σ ≤ 1e−12 rel | R run |
| M-PHSEG-2 | 1000–15 000 SNPs, ≥ 3 Mb gap, < 50-SNP chromosome, options 3 / 1 / 2 | t2 | 6 segments each option | R run |
| M-PHSEG-3 | `runPcfSubset` ≥ 15 000 SNPs | t3 (16 000) | 5 segments | R run |
| M-PHSEG-4 | prior breakpoints, gamma 5, kmin 5; no_segmentation | t4 | 4 segments; 1 segment 0.63947307586669999 | R run |
| M-PHSEG-5 | noisy track (sd 0.1, 7-SNP blocks) | t5 | 4 segments at R's breakpoints | R run |
| M-PHSEG-6 | end-to-end `BuildBattenbergSegments` → `FitSubclonalCopyNumberWithSnpTest` | e1 (t1, ρ 0.7, ψ 2.6, Inf probe), e2 (t2, ρ 0.85, ψ 3) | 9/9 segments: extents, LogR, pval, states, fractions | R run |
| S-PHSEG-7 | NA dropped; < 50 SNPs ⇒ mean phase, median segment BAF | 4 SNPs | phased (0.8, 0.7, 0.1), BAFseg 0.7 | `run_pcf` code |
| S-PHSEG-8 | option 3 keeps mean when median ∈ {0, 1} | (1, 1, 0.7) | 0.9 (option 1: 1) | `run_pcf` code |
| S-PHSEG-9 | `BuildBattenbergSegments`: no finite logR ⇒ 0; chromosome change splits | 3 rows | 2 segments, LogR 0 | `determine_copynumber` code |
| C-PHSEG-1 | invalid arguments | null, BAF 1.2, kmin 0, phasekmin 15, gamma NaN, option 4, BAFseg 1.5 | ArgumentNullException / ArgumentException / ArgumentOutOfRangeException | contract |
| M-STAT-Q1 | `SampleQuantileType7` | 7-value sample, p ∈ {0, 0.05, 0.25, 0.5, 0.85, 0.88, 1}; {3,1,4,1,5} p 0.85 | R `quantile` ≤ 1e−15 | R run |

## 15. FIN-B24 F41 — Battenberg solutions A–F, SDfrac and bootstrap CIs (2026-10-10)

File `OncologyAnalyzer_BattenbergBootstrap_Tests.cs`. Expected values = `set.seed(seed); determine_copynumber(...)`
sourced verbatim (R 4.3.3), every column %.17g (Evidence § F41).

| ID | Test | Input | Expected | Evidence |
|----|------|-------|----------|----------|
| M-BOOT-1 | single segments, all 60 columns | s1, s1n500, s1n20, s2, s3, s4, s5 (seeds 7, 3, −12345, 2³¹−1, 0) | CNs exact, τ/SDfrac/SDfrac_BS/bounds ≤ 1e−12; NA solution F last; NA bounds for noperms < 975 / < 25; 1 SNP ⇒ SDfrac NA, SD_BS 0 | R run |
| M-BOOT-2 | RNG stream across segments | F40 tracks e1 (seed 4711), e2 (seed 99) | every sub-clonal row ≤ 1e−12; clonal rows no solutions | R run |
| M-BOOT-3 | statistical agreement, seeds disjoint from R | s4, 100 C# seeds vs R 400 seeds | means of SD_BS / bounds within 5 MC standard errors; SDfrac seed-independent | R run |
| S-BOOT-4 | reproducibility, clonal ⇒ empty, guards | s1 twice; within-maxdist segment; permutations 0, null, ρ 0 | identical; Solutions empty; ArgumentOutOfRange / ArgumentNull | contract |

## 16. FIN-B24 F53 — multi-sample `ascat.asmultipcf` (2026-10-10)

File `OncologyAnalyzer_AscatAsMultiPcf_Tests.cs` (+ S-ASPCF-G4 in `OncologyAnalyzer_AscatGermlineHomozygous_Tests.cs`).
Expected values = `ascat.asmultipcf(...)` sourced verbatim (R 4.3.3) on the MINSTD cohorts C1–C5 (Evidence § F53).

| ID | Test | Input | Expected | Evidence |
|----|------|-------|----------|----------|
| M-AMPCF-1..14 | `SegmentAlleleSpecificAsMultiPcf_MatchesAscatAsMultiPcf` | C1 (default, no refine, per-track wsample), C2, C3 (default, no refine, wsample, penalty 25, penalty 0.001 ⇒ ladder), C5 (default, wsample c(1, 4)), C4 6500 probes (exact, fast, fast + wsample) | per-sample segment extents / probe / het counts identical; logR and BAF ≤ 1e−12 (observed ≤ 7.8e−16) | R run |
| S-AMPCF-1 | refinement removes unsupported breakpoints | C1 | joint: breaks 400, 900 in both; refined: S1 only 900, S2 only 400 | R run |
| S-AMPCF-2 | R quirk: zero-BAF-weight first block | C1 chr4 (homozygous only) | refine ⇒ 1 level; no refine ⇒ 150 + 150, −1.4828774511765428 | R run |
| S-AMPCF-3 | output shape, default genotypes | C2 | segments tile 1437 probes; BAF NaN exactly at homozygous probes; null genotypes ≡ all heterozygous | contract |
| C-AMPCF-1 | guards | 1 sample, single-probe chromosome (R errors), differing probes, genotype / wsample count, wsample ≤ 0, penalty 0 / NaN, het BAF 1.5, NaN logR | `ArgumentException` / `ArgumentOutOfRangeException` / `ArgumentNullException` | R run / contract |
| S-ASPCF-G4 | F36 last-window `x[n:(n−1)]` | single sample, chr2 of 32 probes | 200 / 100 / 16 / 16 segments, R levels ≤ 1e−14 | R run |

