# Test Specification: ONCO-PURITY-001

**Test Unit ID:** ONCO-PURITY-001
**Area:** Oncology
**Algorithm:** Tumor Purity Estimation
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Antonello et al. (2024), CNAqc, *Genome Biology* 25(1):38 | 1 | https://doi.org/10.1186/s13059-024-03170-5 | 2026-06-14 |
| 2 | CNAqc package vignette (Caravagna lab) | 3 | https://caravagnalab.github.io/CNAqc/articles/CNAqc.html | 2026-06-14 |
| 3 | Carter et al. (2012), ABSOLUTE, *Nat Biotechnol* 30(5):413–421 | 1 | https://doi.org/10.1038/nbt.2203 | 2026-06-14 |
| 4 | Shen & Seshan (2016), FACETS, *NAR* 44(16):e131 | 1 | https://doi.org/10.1093/nar/gkw520 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Expected VAF of a clonal mutation: `v = mπ / [2(1−π) + π(n_A+n_B)]`, where m = multiplicity, π = purity, n_A+n_B = total copy number — CNAqc vignette (source 2), CNAqc paper (source 1).
2. Copy-neutral diploid heterozygous SNV (m=1, n_tot=2): v = π/2 ⇒ purity ρ = 2·VAF — derived from key point 1; confirmed by CNAqc 60%/30% example (source 1).
3. Inverting the general relation for purity: π = 2v / [m − v(n_tot − 2)] — algebraic inversion of key point 1 with c = 1 (see §2.2).
4. 2:1 segment (n_tot=3) at π=1: clonal VAF peaks at 1/3 (m=1) and 2/3 (m=2) — CNAqc paper (source 1).
5. Normal cells contribute exactly 2 copies weighted by (1−π); FACETS mixing `m* = mΦ+(1−Φ)` independently confirms the 2(1−π) term (source 4).
6. Invariant 0 ≤ purity ≤ 1 — checklist invariant; intrinsic to the fraction-of-tumor-cells definition (source 3).

### 1.3 Documented Corner Cases

- Multiplicity ambiguity on amplified segments (n_tot>2): VAF alone is insufficient; copy-neutral diploid loci are the robust case — CNAqc (source 1).
- Subclonal mutations (c<1) depress VAF and underestimate purity if treated as clonal — CNAqc (source 1).
- Purity < 0.1 below detection; no informative heterozygous variants; high stromal contamination — checklist by-area definition.

### 1.4 Known Failure Modes / Pitfalls

1. Using a subclonal or amplified-segment VAF in the diploid ρ=2·VAF formula overestimates/underestimates purity — CNAqc (source 1).
2. VAF > 0.5 under the diploid heterozygous model implies purity > 1 (impossible) — indicates wrong copy state or LOH; must be rejected — derived from §2.2 domain.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `EstimatePurityFromVAF(IEnumerable<VariantObservation>)` | OncologyAnalyzer | Canonical | Copy-neutral diploid heterozygous model: per-variant ρ = 2·VAF, aggregated by median. |
| `EstimatePurity(IEnumerable<PurityVariant>)` | OncologyAnalyzer | Canonical | Allele-specific model: inverts v = mπ/[2(1−π)+π·n_tot]; aggregated by median. |
| `EstimatePurityFromVaf(double vaf)` | OncologyAnalyzer | Delegate | Single-VAF closed form ρ = 2·VAF; smoke only. |
| `AnalyzePurityPeaks(IEnumerable<PurityPeakMutation>, double, PurityPeakOptions?)` | OncologyAnalyzer | Canonical | CNAqc 1.1.5 `analyze_peaks_common` peak-based purity QC (FIN-B24 F34). |
| `StatisticsHelper.BandwidthNrd0` / `GaussianKernelDensity` / `PeakPick` | StatisticsHelper (Infrastructure) | Canonical (shared) | R `bw.nrd0` / `density.default` / `peakPick::peakpick` (FIN-B24 F33). |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ purity ≤ 1 | Yes | Checklist invariant; fraction-of-cells definition (source 3) |
| INV-2 | For m=1, n_tot=2: purity = 2·VAF exactly | Yes | §1.2 pt 2; CNAqc 60%/30% example (source 1) |
| INV-3 | Inversion of the expected-VAF relation recovers the purity that generated the VAF | Yes | §1.2 pt 3 (algebraic inverse of source 1/2 formula) |
| INV-4 | A VAF of 0 yields purity 0; the estimator is monotone non-decreasing in VAF for fixed m, n_tot | Yes | ρ=2v closed form (source 1) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | VAF-only diploid het, VAF 0.30 | Single clonal het SNV at VAF 0.30 | purity = 0.60 | CNAqc 60%/30% example (source 1) |
| M2 | VAF-only boundary VAF 0.50 | het SNV at VAF 0.50 | purity = 1.0 | ρ=2·VAF (source 1) |
| M3 | VAF-only VAF 0.00 | absent variant | purity = 0.0 | ρ=2·VAF; INV-4 |
| M4 | VAF-only median of several variants | VAFs {0.10,0.15,0.30} (purities {0.20,0.30,0.60}) | purity = 0.30 (median) | §1.2 pt 2 + median aggregation |
| M5 | `EstimatePurityFromVaf` single delegate | vaf 0.275 | purity = 0.55 | ρ=2·VAF; CNAqc 0.275↔0.55 band (source 1) |
| M6 | Allele-specific diploid het VAF 0.30 | PurityVariant(vaf 0.30, m 1, n_tot 2) | purity = 0.60 (agrees with M1) | CNAqc example (source 1) |
| M7 | Allele-specific 2:1 m=1 VAF 1/3 | PurityVariant(vaf 1/3, m 1, n_tot 3) | purity = 1.0 | CNAqc 2:1 peak 33% (source 1) |
| M8 | Allele-specific 2:1 m=2 VAF 2/3 | PurityVariant(vaf 2/3, m 2, n_tot 3) | purity = 1.0 | CNAqc 2:1 peak 66% (source 1) |
| M9 | Allele-specific general non-trivial | PurityVariant(vaf 1/6, m 1, n_tot 4): π = 2·(1/6)/[1 + (1/6)(2−4)] = (1/3)/(2/3) | purity = 0.5 | inversion §1.2 pt 3 (forward: π=0.5,m=1,n=4 ⇒ v=0.5/3=1/6) |
| M10 | VAF-only invalid VAF > 0.5 | het SNV VAF 0.6 ⇒ purity>1 | throws ArgumentOutOfRangeException | INV-1 domain |
| M11 | VAF-only empty input | no variants | throws ArgumentException | undefined purity (source 1 corner case) |
| M12 | Allele-specific invalid copy number | PurityVariant n_tot 0 | throws ArgumentOutOfRangeException | formula domain (denominator) |
| M13 | VAF-only invalid VAF range | VAF −0.1 / 1.1 | throws ArgumentOutOfRangeException | VAF ∈ [0,1] |
| M14 | null input | null enumerable | throws ArgumentNullException | standard guard |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Low purity below detection | het SNV VAF 0.02 | purity = 0.04 (no error) | checklist: purity<0.1 |
| S2 | Allele-specific median across variants | mixed PurityVariants | median purity | robustness |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | same input twice | identical result | order-independent / deterministic |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests existed for `EstimatePurity` / `EstimatePurityFromVAF` / `EstimatePurityFromVaf` (the methods are new in this unit). Sibling VAF tests live in `OncologyAnalyzer_CalculateVAF_Tests.cs` (ONCO-VAF-001) and are unrelated.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14, S1, S2, C1 | ❌ Missing | New unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_EstimatePurity_Tests.cs` — all cases for this unit.
- **Remove:** (none)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_EstimatePurity_Tests.cs` | Canonical (all cases) | 21 (incl. TestCase rows) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented (VAF-only + read-count) | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented (median) | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented (delegate) | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented (VAF-only + allele-specific) | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented (both overloads) | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented (both overloads) | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented (both overloads) | ✅ Done |
| 15 | S1 | ❌ Missing | Implemented | ✅ Done |
| 16 | S2 | ❌ Missing | Implemented | ✅ Done |
| 17 | C1 | ❌ Missing | Implemented (determinism) | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `..._DiploidHetVaf30Percent...`, `..._SingleDiploidHet...` |
| M2 | ✅ Covered | `..._ClosedForm_ReturnsTwiceVaf` (0.50) |
| M3 | ✅ Covered | `..._ClosedForm_ReturnsTwiceVaf` (0.00) |
| M4 | ✅ Covered | `..._MultipleVariants_ReturnsMedianPurity` |
| M5 | ✅ Covered | `..._ClosedForm_ReturnsTwiceVaf` (0.275) |
| M6 | ✅ Covered | `..._DiploidHet_AgreesWithVafOnly` |
| M7 | ✅ Covered | `..._TwoToOneSegment_PurePeaks_ReturnOne` |
| M8 | ✅ Covered | `..._TwoToOneSegment_PurePeaks_ReturnOne` |
| M9 | ✅ Covered | `..._GeneralAlleleSpecific_RecoversPurity` |
| M10 | ✅ Covered | `..._VafAbovePointFive_Throws`, `..._CombinationImplyingPurityAboveOne_Throws` |
| M11 | ✅ Covered | `..._Empty_Throws`, `..._EmptyAndNull_Throw` |
| M12 | ✅ Covered | `..._InvalidCopyNumberOrMultiplicity_Throws` |
| M13 | ✅ Covered | `..._VafOutOfRange_Throws` (both overloads) |
| M14 | ✅ Covered | `..._Null_Throws`, `..._EmptyAndNull_Throw` |
| S1 | ✅ Covered | `..._LowPurity_ReturnsSmallPurityNoError` |
| S2 | ✅ Covered | `..._MultipleVariants_ReturnsMedian` |
| C1 | ✅ Covered | `..._RepeatedCalls_AreDeterministic` |

---

### 5.7 Review 2026-09 (B24) additions

| Test | Evidence |
|------|----------|
| `EstimatePurity_CnaqcExpectedVafFun_RoundTripRecoversPurity` (6 karyotypes) | CNAqc `R/equations.R` `expected_vaf_fun`, Python port values |
| `EstimateCcf_RawCcf_EqualsCanonicalPurityCorrectionOverMultiplicity` | dedup lock: EstimateCcf/DeriveMultiplicity route through `AdjustVAFForPurity` |

| `EstimatePurity_MultiplicityAboveTotalCopyNumber_Throws` (v 0.5, m 3, n_tot 2) | CNAqc `expectations_generalised` m ∈ 1..Major (R source); FIN-B24 F24 |
| `EstimatePurity_ExactClonalPeakAtFullPurity_ReturnsExactlyOne` ((0.2,1,5), (0.4,2,5) → 1.0) | CNAqc `expected_vaf_fun(1,4,1,1)` = 0.20000000000000001, `(2,3,2,1)` = 0.40000000000000002 in R; raw inversion 1.0000000000000002; FIN-B24 F25 |
| `EstimatePurity_SlightlyAboveFullPurityPeak_StillThrows` (v 0.2·(1+1e-9), m 1, n_tot 5 → π 1.0000000025000004) | tolerance guard; FIN-B24 F25 |

Property `KnownPurityVariantGen` (Properties/OncologyProperties.cs) now draws π ∈ (0, 1] (cap 950‰ → 1000‰) since the π = 1 boundary no longer throws (FIN-B24 F25; property class passes, 422 tests).

Fixed (FIN-B24 F24/F25, code in B22-owned `OncologyAnalyzer.SomaticCalling.cs` under orchestrator special permission).

### 5.8 FIN-B24 F33/F34 — CNAqc peak-based purity QC

Data: `TestData/CNAqc/cnaqc_D{1,3,4}.tsv` — seeded R simulations (D1 π 0.7: 1:1 ×400, 2:1 ×250, 2:2 ×150, 1:0 ×120, 3:1 ×60; D3 π 0.45: 1:1 ×500 clonal + ×200 CCF 0.35, 2:0 ×150, + 60 true-4:0 m = 4 mutations labelled 2:0; D4 π 0.3: 1:1 ×150, 2:1 ×90, 1:0 ×110). Reference: CNAqc 1.1.5 (4b7cea4a) `analyze_peaks_common` sourced in R 4.3.3 (prefix-strip + `easypar::run` shims), density lattice R ≥ 4.4 unless "old".

| ID | Test | Expected (CNAqc R) |
|----|------|--------------------|
| P1 | D1 π 0.7 default | λ 0.0023756289876209163, PASS (all karyotypes PASS); full match + peak tables |
| P2 | D1 π 0.7 `LegacyDensityCoordinates` (R 4.3.3 native `density`) | λ identical, peak y 8.64 vs 8.63 (1:0) |
| P3 | D1 π 0.5 (wrong purity) | λ −0.29645326597409133, FAIL (all karyotypes FAIL) |
| P4 | D1 π 0.7 `KernelAdjust` 0.5 | λ 0.023857243903529196, PASS; 1:1 6 KDE peaks incl. one discarded (y 0.31 ≤ 8.19/20) |
| P5 | D1 π 0.7 ε 0.01, tol 0.005 | 1:0 FAIL, sample PASS (weight 0.87 vs 0.13) |
| P6 | D1 π 0.62 | λ −0.11362481021020236 (|λ| > ε) yet PASS; 2:1 m = 2 unmatched but 2:1 PASS (lead row = highest counts_per_bin) |
| P7 | D3 π 0.45 / `MinVaf` 0.12 | λ 0.014615384615384629 / 0.030000000000000027 (subclonal tail removed changes weights) |
| P8 | D4 π 0.3: default / size 0.35 / min count 80 / min count 200 | λ 0.0039996443291489122 / 0 (only 1:1) / −0.010660215188956853 (2:1 included) / not analysed (`Pass` null, λ NaN) |
| P9 | Rightmost (legacy `peak_detector` matching) on D1, D3 | same peaks as Closest; rows in descending expected-peak order |
| P10 | Mixture peaks = R BMix `B.params` (set.seed 7): D1, D3 Closest, D3 Rightmost | λ 0.0045358591466179345, 0.0037889790909820253, −0.13824761798675933; D3 2:0 Closest PASS vs Rightmost FAIL (0.45 → BMix 0.618 peak) |
| P11 | Σ karyotype scores = λ; Σ weights = 1; counts 120/400/250/150 | identity |
| P12 | MinKaryotypeSize denominator includes non-simple karyotypes (1:0 = 120/980) | kept at 120/980, dropped at +1e−9 |
| P13 | invalid arguments (null, π ∉ (0,1], ε 0, size 1, adjust 0, min count −1, karyotype 3:1, VAF 1.2, Major −1) | throws |
| K1 | `BandwidthNrd0` ×3 datasets + 4 fall-back cases | R `bw.nrd0` |
| K2 | `GaussianKernelDensity` ×6 (dataset, adjust, coords) × 9 grid points | R `density` (≤ 1e−12 rel.) |
| K3 | `PeakPick` on R densities (neighlim 1–5) and synthetic series (neighlim 0/1/8/12/20; deriv.lim, peak.min.sd, npos) | peakPick 0.11 |

Tests: `OncologyAnalyzer_AnalyzePurityPeaks_Tests` (21), `StatisticsHelper_GaussianKde_Tests` (26).

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | VAF-only estimator assumes clonal, heterozygous (m=1), copy-neutral diploid (n_tot=2) SNVs | `EstimatePurityFromVAF`, `EstimatePurityFromVaf` |
| 2 | Multiple per-variant estimates aggregated by median (robust central estimator) | `EstimatePurityFromVAF`, `EstimatePurity` |

---

## 7. Open Questions / Decisions

1. Decision: the VAF-only estimator deliberately fixes the copy-number model to copy-neutral diploid heterozygous (the textbook ρ=2·VAF case); the allele-specific `EstimatePurity` covers other copy states via explicit (m, n_tot) inputs. Both are fully source-derived; no correctness-affecting assumption remains.
