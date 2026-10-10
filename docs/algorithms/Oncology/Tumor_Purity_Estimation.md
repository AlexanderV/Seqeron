# Tumor Purity Estimation

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-PURITY-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-10-10 |

## 1. Overview

Tumor purity ρ (also π) is the fraction of cells in a bulk sequencing sample that are tumour cells, the remainder being normal/stromal cells. These estimators recover ρ from the variant allele frequencies (VAFs) of clonal somatic mutations by inverting the closed-form expected-VAF relation that links VAF to purity, mutation multiplicity, and local copy number [1][2]. The computation is exact (deterministic, closed form) for a given (VAF, multiplicity, copy-number) state; the canonical special case — a clonal heterozygous somatic SNV at a copy-neutral diploid locus — gives ρ = 2·VAF [1]. It should be used when somatic SNV calls (and, for non-diploid loci, allele-specific copy-number state) are available; it is not a copy-ratio or B-allele-frequency segmentation method. CNAqc's own purity procedure — the peak-based QC of `analyze_peaks` for simple clonal karyotypes — is ported as `AnalyzePurityPeaks` (§4.4): it does not estimate purity but scores a supplied purity against the VAF peaks it implies [1]. The two companion analyses `analyze_peaks` runs on the same purity — complex clonal karyotypes (`analyze_peaks_general`) and subclonal segments (`analyze_peaks_subclonal`) — are ported as `AnalyzeComplexKaryotypePeaks` and `AnalyzeSubclonalPurityPeaks` (§4.5).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A bulk tumour sample mixes tumour cells (fraction ρ) with normal diploid cells (fraction 1−ρ). For a somatic mutation present in m copies of the tumour genome on a segment of total copy number n_tot, the alternate-allele reads come only from tumour cells, while the total reads come from both populations. The expected VAF therefore depends jointly on ρ, m, and n_tot [1][2]. ABSOLUTE [3] and FACETS [4] use the same mixture to convert allelic fractions to per-cancer-cell allele counts and to estimate purity/ploidy.

### 2.2 Core Model

CNAqc gives the expected VAF of a clonal mutation (cancer cell fraction c = 1) as [1][2]:

> v = m·π / [ 2(1−π) + π·(n_A + n_B) ]

where m is the mutation multiplicity, π the tumour purity, and n_A + n_B = n_tot the tumour total copy number; the term 2(1−π) is the contribution of the healthy diploid normal cells and π·n_tot is the tumour contribution [2]. FACETS independently confirms the normal-diploid term: it mixes a normal (1,1) genotype with the aberrant (m,p) genotype at cellular fraction Φ via m* = mΦ + (1−Φ) [4].

Inverting for purity:

> π = 2·v / [ m + v·(2 − n_tot) ]

For a clonal **heterozygous** SNV at a **copy-neutral diploid** locus (m = 1, n_tot = 2) this collapses to v = π/2, i.e. **ρ = 2·v** [1]. CNAqc's worked example states a real purity of 60% corresponds to VAF 30% (and the 55–65% purity band to the 27.5–32.5% VAF band) [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Mutations are clonal (cancer cell fraction c = 1) | A subclonal mutation (c<1) has depressed VAF (v ∝ c); purity is underestimated [1] |
| ASM-02 | For the VAF-only estimator: m = 1 and n_tot = 2 (copy-neutral diploid heterozygous) | On amplified/LOH segments ρ = 2·VAF is wrong; use the allele-specific `EstimatePurity` with explicit (m, n_tot) [1] |
| ASM-03 | Normal cells are diploid (2 copies) | The 2(1−π) term and the closed form no longer hold |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ purity ≤ 1 | Purity is a fraction of cells [3]; inputs yielding ρ outside [0,1] are rejected |
| INV-02 | For m=1, n_tot=2: purity = 2·VAF exactly | v = π/2 special case of §2.2 [1] |
| INV-03 | Inversion recovers the purity that generated the VAF | π = 2v/[m+v(2−n_tot)] is the exact algebraic inverse of the §2.2 relation |
| INV-04 | VAF = 0 ⇒ purity = 0; estimate is monotone non-decreasing in VAF for fixed m, n_tot | ρ = 2v closed form [1] |

### 2.5 Comparison with Related Methods

| Aspect | This estimator (VAF/closed-form) | ABSOLUTE / FACETS |
|--------|----------------------------------|-------------------|
| Input | Clonal SNV VAFs (+ optional allele-specific CN) | Genome-wide segmented copy ratios + BAF + SNVs |
| Output | Point purity from the closed-form relation | Joint purity + ploidy + absolute CN by model fitting |
| Multiplicity | Supplied per variant | Inferred jointly |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| variants (`EstimatePurityFromVAF`) | `IEnumerable<VariantObservation>` | required | Clonal heterozygous diploid somatic SNVs | non-null, non-empty; valid read counts; each VAF ≤ 0.5 |
| vaf (`EstimatePurityFromVaf`) | `double` | required | Single clonal het diploid SNV VAF | ∈ [0, 0.5] |
| variants (`EstimatePurity`) | `IEnumerable<PurityVariant>` | required | Clonal SNVs with VAF, multiplicity m, total CN n_tot | non-null, non-empty; 1 ≤ m ≤ n_tot; VAF ∈ [0,1] |
| mutations (`AnalyzePurityPeaks`) | `IEnumerable<PurityPeakMutation>` | required | Somatic mutations with VAF and segment karyotype Major:minor (all karyotypes) | non-null; VAF ∈ [0,1]; Major, minor ≥ 0 |
| purity (`AnalyzePurityPeaks`) | `double` | required | Purity π to QC | ∈ (0, 1] |
| options (`AnalyzePurityPeaks`) | `PurityPeakOptions?` | CNAqc defaults | karyotypes {1:0,1:1,2:0,2:1,2:2}, `MinKaryotypeSize` 0, `MinAbsoluteKaryotypeMutations` 100, `PurityError` ε 0.05, `VafTolerance` 0.015, `KernelAdjust` 1, `MatchingStrategy` Closest, `MinVaf` 0, `MixturePeaks` null, `LegacyDensityCoordinates` false, `FitMixturePeaks` false (CNAqc itself always fits BMix; needs `AlternateReads`/`Depth` on the mutations), `BootstrapCount` 1 (`n_bootstrap`), `Seed` 0 (R `set.seed`) | ε ∈ (0,1); size ∈ [0,1); adjust > 0; karyotypes ⊆ simple set; bootstrap ≥ 1; not both `FitMixturePeaks` and `MixturePeaks` |
| successes, trials, seed (`FitBinomialMixture`) | `IReadOnlyList<int>` ×2, `int` | required | NV, DP per mutation; R `set.seed` value; optional component counts (default 1:4) | 1 ≤ DP, 0 ≤ NV ≤ DP; K ≥ 1 |
| mutations, purity, options (`AnalyzeComplexKaryotypePeaks`) | as `AnalyzePurityPeaks` | CNAqc defaults | uses `Karyotypes` (gate only), `MinAbsoluteKaryotypeMutations` (`n_min`), `PurityError` (`epsilon`), `KernelAdjust`, `MinVaf`, `LegacyDensityCoordinates` | as above |
| segments (`AnalyzeSubclonalPurityPeaks`) | `IEnumerable<SubclonalPeakSegment>` | required | subclonal segments: karyotype 1 at `Ccf`, karyotype 2 at 1 − `Ccf`, segment VAFs | both karyotypes simple; CCF ∈ (0,1); VAF ∈ [0,1] |
| options (`AnalyzeSubclonalPurityPeaks`) | `SubclonalPeakOptions?` | `analyze_peaks` values | `MinMutations` 100 (n > it), `Epsilon` 0.05, `KernelAdjust` 1, `StartingState` 1:1, `Seed` 0 (R `set.seed` for the mutation identifiers), `LegacyDensityCoordinates` false | ε ∈ (0,1); start Major ≥ minor ≥ 0, Major ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | Estimated tumour purity ρ ∈ [0, 1]; median of per-variant estimates for the collection overloads |
| (return, `AnalyzeComplexKaryotypePeaks`) | `ComplexKaryotypePeakAnalysis` | `Ran` (gate), per-karyotype `ComplexKaryotypePeakResult` in CNAqc summary order (matched / mismatched expected peaks, `MatchedProportion`, `Pass` = prop ≥ 0.5, KDE, data peaks) and all `ComplexKaryotypeExpectedPeak`s (m, expected VAF, matched) |
| (return, `AnalyzeSubclonalPurityPeaks`) | `IReadOnlyList<SubclonalSegmentPeakResult>` | per analysed segment: KDE, data peaks, every model's expected peaks (`SubclonalExpectedPeak`: model id, linear/branching, copies in each clone, genotypes, shared/private, expected VAF, matched), `Rankings` and `BestModels` (CNAqc decision table) |
| (return, `AnalyzePurityPeaks`) | `PurityPeakAnalysis` | `Score` λ = Σ weight·offset (purity units; CNAqc prints it as "Purity correction"), `Pass` (sample QC; null = no karyotype passed the filters), per-karyotype `PurityPeakKaryotype` (n, weight, score, PASS/FAIL, KDE, data peaks) and per-expected-peak `PurityPeakMatch` (expected peak, δ, matched peak, offset_VAF, offset, weight, matched) |

### 3.3 Preconditions and Validation

Null collections throw `ArgumentNullException`; empty collections throw `ArgumentException` (purity undefined). A VAF outside [0,1] throws `ArgumentOutOfRangeException`; for the diploid model a VAF > 0.5 (implying ρ > 1) throws `ArgumentOutOfRangeException`. For the allele-specific overload, m < 1, n_tot < 1, m > n_tot, or any (VAF, m, n_tot) combination yielding ρ outside [0,1] (including a non-positive denominator), throws `ArgumentOutOfRangeException`; a computed ρ in (1, 1 + (n_tot+4)·ε] is IEEE rounding of an exact π = 1 peak and is clamped to 1.0. Read counts are validated via `CalculateVAF` (alt/total) as in ONCO-VAF-001. `AnalyzePurityPeaks`: null mutations → `ArgumentNullException`; π ∉ (0, 1], ε ∉ (0, 1), size ∉ [0, 1), negative minimum count, adjust ≤ 0, negative/infinite tolerance, NaN `MinVaf`, a VAF ∉ [0, 1] or a negative allele copy number → `ArgumentOutOfRangeException`; a requested karyotype outside {1:0, 1:1, 2:0, 2:1, 2:2} → `ArgumentException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. For each variant compute its VAF (read-count overload) or read its supplied VAF.
2. Map the VAF to a per-variant purity: ρ = 2·v (diploid model) or ρ = 2v/[m + v(2 − n_tot)] (allele-specific).
3. Validate ρ ∈ [0, 1].
4. Aggregate per-variant purities by their median.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Normal total copy number fixed at 2 (autosomal diploid normal) [2][4].
- Heterozygous-diploid factor: ρ = 2·VAF (m=1, n_tot=2) [1].
- Median chosen over mean for robustness to subclonal/outlier VAFs (aggregation policy; does not alter the single-variant formula).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| EstimatePurity / EstimatePurityFromVAF | O(n log n) | O(n) | n = #variants; dominated by the median sort. O(1) per variant. |
| EstimatePurityFromVaf | O(1) | O(1) | single closed-form evaluation |
| AnalyzePurityPeaks | O(n + K·512²) | O(n) | per karyotype: one 512-point KDE (direct convolution) + peakPick ×5 |
| AnalyzeComplexKaryotypePeaks | O(n + K·512²) | O(n) | as above, no mixture peaks |
| FitBinomialMixture | O(8·(100·I·n·K + n·K)) | O(n·K) | 8 grid runs × 100 k-means starts (Hartigan–Wong, I ≤ 10 iterations) + one E/M step; long-double sums emulated with big integers |
| AnalyzeSubclonalPurityPeaks | O(S·(n_s + 512² + E)) | O(n + E) | E = evolution states explored (breadth-first, ploidy ≤ 2·target ploidy; a few hundred for simple karyotypes) |

### 4.4 CNAqc peak-based purity QC (`AnalyzePurityPeaks`)

Port of CNAqc 1.1.5 `analyze_peaks` → `analyze_peaks_common` (caravagnalab/CNAqc `R/peak_algorithms.R`, `R/vaf_functions.R`, `R/equations.R`, `R/peak_detector.R`) [1]:

1. Keep mutations with VAF > `min_VAF`; count n_K per karyotype over **all** karyotypes (N = Σ n_K). Analyse K ∈ `karyotypes` with n_K ≥ `min_absolute_karyotype_mutations` and n_K/N ≥ `min_karyotype_size` (mutations of all segments with the same karyotype are pooled). None → no analysis.
2. Data peaks per K: Gaussian KDE of the VAFs (R `density`, `bw.nrd0`·`kernel_adjust`, n = 512, cut = 3; `StatisticsHelper.GaussianKernelDensity`); union of `peakPick::peakpick(neighlim = 1..5)` maxima (`StatisticsHelper.PeakPick`); x, y rounded to 2 decimals (R `round`), distinct x, kept in [0, 1]; `counts_per_bin` = `hist(VAF, breaks = seq(0, 1, 0.01))` count at bin round(100·x); discarded ⇔ y ≤ max(y)/20. CNAqc also appends BMix Binomial-mixture means (`bmixfit(K.Binomials = 1:4)`, §4.6) snapped to the KDE grid — fitted with `FitMixturePeaks` or supplied via `MixturePeaks`; `n_bootstrap` > 1 adds resampled peaks (§4.6).
3. Expected peaks for m ∈ {1, Major}: v_m = m·π / (2(1−π) + π·(Major+minor)) (`ascat()` in `vaf_functions.R`); band δ_m = 2·m·ε / (2 + π·(ploidy−2))² (`delta_vaf_karyo`).
4. Matching (`Closest`, the only strategy 1.1.5 runs): each v_m takes the nearest non-discarded data peak x. (`Rightmost`, legacy `peak_detector`: expected peaks descending ↔ highest-x data peaks descending, padded with the rightmost.)
5. offset = 2·m·(v_m − x) / (m + x·(2 − ploidy))² (`compute_delta_purity`, purity units); weight = n_K / Σ n_analysed; matched ⇔ max(x − tol, v_m − δ_m) ≤ min(x + tol, v_m + δ_m) (`overlap_bands`).
6. Karyotype QC = `matched` of its row with the largest `counts_per_bin` (stable, first on ties); sample QC = PASS/FAIL class with the larger Σ weight over rows (tie → FAIL); λ = Σ weight·offset.

CNAqc does **not** compare λ with ε (ε only sizes the bands), the `p_binsize_peaks` argument is unused in 1.1.5, and no corrected purity is proposed — `print.cnaqc` shows λ as "Purity correction". Example (R-locked, seeded dataset D1, true π = 0.7): at π = 0.7 λ = 0.00238, PASS; at π = 0.5 λ = −0.296, FAIL; at π = 0.62 the 1:0 karyotype fails but the sample passes by weight (0.87) although |λ| = 0.114 > ε.

### 4.5 Complex and subclonal karyotypes (`AnalyzeComplexKaryotypePeaks`, `AnalyzeSubclonalPurityPeaks`)

Ports of CNAqc 1.1.5 `analyze_peaks_general` and `analyze_peaks_subclonal` (`R/peak_algorithms.R`) with `expectations_generalised` / `expectations_subclonal` (`R/equations.R`), called by `analyze_peaks` after the simple-karyotype QC with `n_min = min_absolute_karyotype_mutations`, `epsilon = purity_error` [1]:

- **Complex clonal karyotypes.** `analyze_peaks` runs the step only if some karyotype outside `karyotypes` has n > `min_absolute_karyotype_mutations` (strict); it then analyses every karyotype other than the fixed simple five with n ≥ `n_min` (inclusive), in first-appearance order. Data peaks = the simple-karyotype KDE detector (no BMix). Expected peaks for m = 1..max(Major, minor, 1): m·π / (2(1−π) + π·ploidy). An expected peak is matched iff some data peak — discarded ones included — satisfies |x − v_m| < ε. Summary: matched / mismatched counts, prop = matched / total, ordered by descending prop (ties by the "Major:minor" string); `analyze_peaks` labels the karyotype's segments and mutations QC PASS iff prop ≥ 0.5. No score or sample verdict.
- **Subclonal segments** (`cluster_subclonal_CCF = FALSE`; each segment with n > `n_min` on its own). Evolution models from the starting state (1:1 = alleles A1 B1 by default): breadth-first single-allele amplification / deletion and whole-genome doubling, ploidy capped at 2·(target ploidy), until the target karyotype appears; distinct allele sets kept; every reached state gives each allele a new mutation. Models: branching (both clones from the start) and linear in both directions (clone 1 → clone 2 only if clone 1 has no LOH or both have LOH; symmetric for clone 2 → clone 1 with cell fraction 1 − CCF). For a clone pair each mutation with n1, n2 copies gives the peak π·(n1·CCF + n2·(1 − CCF)) / (2(1−π) + π·(CCF·ploidy1 + (1 − CCF)·ploidy2)); distinct values, ascending; clone pairs with identical peak vectors (R `paste`, 15 significant digits) are dropped. Matched ⇔ |x − peak| ≤ ε for some KDE data peak; each model scores matched / #peaks; the top-scoring models are the decision (CNAqc `summary`). CNAqc names mutations by `sample(LETTERS, 8, replace = TRUE)`; the port draws them from the R RNG (`Seed` = `set.seed`), so identifiers — and which of several equal-VAF mutations is reported — are reproduced; peaks and decisions do not depend on the seed.
- CNAqc quirks kept: a LOH starting state cannot reach a karyotype without LOH → `cli_abort` for that segment, caught by `easypar` → the segment keeps its data peaks but has no model; `evolve()` loops forever when every move exceeds the ploidy cap (e.g. start 2:2 → 1:0) — the port throws `InvalidOperationException` instead.

Example (R-locked, seeded segment S1 chr1, 2:1 at CCF 0.6 / 1:1, π = 0.7): three models; branching A1B1 → A1A2B1 | A1B1 and linear A1B1 → A1A2B1 → A2B1 both match 3 of 4 peaks (0.116, 0.174, 0.289, 0.463; data peaks 0.14, 0.45) and are reported; linear A1B1 → A1B1 → A1A2B1 scores 2/3.

### 4.6 BMix mixture peaks and `n_bootstrap` (`FitBinomialMixture`, `FitMixturePeaks`, `BootstrapCount`)

CNAqc's default `analyze_peaks` runs, for every simple karyotype, `combined_peak_detector` = KDE peaks + `mixture_peak_detector` → `BMix::bmixfit(data.frame(NV, DP), K.Binomials = 1:4, K.BetaBinomials = 0)` (caravagnalab/BMix `R/bmixfit.R`, `R/bmixfit_EM.R`). Complex and subclonal analyses use KDE peaks only. Port (`FitBinomialMixture`, reproducing R's random stream after `set.seed`):

1. Grid K = 1..4 × `samples = 2` runs (K ascending); each run (`runner`, up to 14 attempts on error): `kmeans(NV/DP, K, nstart = 100)` — R `stats::kmeans` with Hartigan–Wong (AS 136, `kmns.f`, `iter.max = 10`; MacQueen when K = 1), every start = `sample.int` of K distinct values, the start with the smallest Σ withinss kept (strict <); K larger than the number of distinct frequencies is an R error and that K drops out.
2. `sample()` permutes the K centres; means = centres + `runif(K, −0.025, 0.025)`, redrawn until all lie in (0, 1); mixing proportions = k-means cluster shares.
3. EM: BMix's loop stops after its **first** iteration — its test `NLL_new − NLL_old < epsilon` starts from `NLL_old = .Machine$integer.max`, so the change is always negative. One E-step: z_ik ∝ π_k·Binom(NV_i | DP_i, p_k) (R `dbinom` log, log-sum-exp), NLL = −Σ_i log Σ_k; one M-step: p_k = Σ z·NV / Σ z·DP, π_k = Σ z / N.
4. BIC = 2·NLL + ln N·2K, ICL = BIC − Σ z·ln z (NaN terms dropped); `which.min(ICL)` over the grid. CNAqc snaps each ICL-best mean to the nearest KDE grid x (never discarded; counts from the same histogram).
5. `n_bootstrap` = B > 1 (`simple_peak_detector`): B resamples `sample_n(replace = TRUE)` of the karyotype's VAFs, their KDE peaks pooled, `distinct(x)`, heights re-read from the full-data density (`phase_to_density`), new x values appended after the full-data peaks. With BMix (`mixture_peak_detector`): the full-data fit still runs (its random numbers are drawn) but its peaks are replaced by those of B resampled fits (each snapped to its own resample's density), `distinct(x)`, heights from the density of one more resample. Stream order per karyotype: KDE resamples, BMix full fit, BMix resamples, the extra resample; karyotypes in "Major:minor" order. `AnalyzeComplexKaryotypePeaks` / `AnalyzeSubclonalPurityPeaks` bootstrap the KDE detector the same way (subclonal: after all model identifiers).

R sums run in 80-bit long double (`sum`, `colSums`) — emulated exactly by `StatisticsHelper.ExtendedPrecisionSum`; `dbinom` = `StatisticsHelper.BinomialLogDensity`. R's `%*%` uses the system BLAS (OpenBLAS here); the port uses R's `matprod = "internal"` long-double dot product, bit-identical to R with that option and within 1e−15 relative of OpenBLAS. Example (R-locked): D1 2:1 (n = 250, set.seed 7) → ICL 3481.4, 3503.0, **1901.15**, 1911.4, 2001.6, 1999.0, 2075.3, 2118.0 → K = 2, means 0.51668 / 0.26290 (expected 0.519 / 0.259 at π 0.7).

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.SomaticCalling.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.SomaticCalling.cs) (record `PurityVariant` and constants in `OncologyAnalyzer.cs`)

- `OncologyAnalyzer.EstimatePurityFromVAF(IEnumerable<VariantObservation>)`: median of ρ = 2·VAF over clonal het diploid SNVs.
- `OncologyAnalyzer.EstimatePurityFromVaf(double)`: single-VAF closed form ρ = 2·VAF.
- `OncologyAnalyzer.EstimatePurity(IEnumerable<PurityVariant>)`: median of the allele-specific inversion ρ = 2v/[m+v(2−n_tot)].
- `OncologyAnalyzer.AnalyzePurityPeaks(IEnumerable<PurityPeakMutation>, double purity, PurityPeakOptions? options = null)` ([OncologyAnalyzer.CopyNumberPloidy.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.CopyNumberPloidy.cs)): CNAqc peak-based purity QC (§4.4); KDE / peak picking via `StatisticsHelper.GaussianKernelDensity` / `PeakPick` (Infrastructure).
- `OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(IEnumerable<PurityPeakMutation>, double purity, PurityPeakOptions? options = null)` (same file): CNAqc `analyze_peaks_general` (§4.5).
- `OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(IEnumerable<SubclonalPeakSegment>, double purity, SubclonalPeakOptions? options = null)` (same file): CNAqc `analyze_peaks_subclonal` + `expectations_subclonal` (§4.5); mutation identifiers from the R Mersenne-Twister port (`RMersenneTwister`, F41).
- `OncologyAnalyzer.FitBinomialMixture(IReadOnlyList<int> successes, IReadOnlyList<int> trials, int seed, IReadOnlyList<int>? componentCounts = null)` (same file): BMix `bmixfit` (§4.6) with private ports of R `kmeans` (Hartigan–Wong `kmns.f`, MacQueen); `StatisticsHelper.ExtendedPrecisionSum` (R long-double `sum`) and `StatisticsHelper.BinomialLogDensity` (R `dbinom(log = TRUE)`) in Infrastructure.

### 5.2 Current Behavior

Collection overloads aggregate per-variant purities by median (lower-mid average for even counts). The VAF-only overload fixes the copy state to copy-neutral diploid heterozygous; non-diploid states are handled only by the allele-specific overload with explicit (m, n_tot). No search/matching is involved, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ρ = 2·VAF for a clonal heterozygous copy-neutral diploid SNV (m=1, n_tot=2) [1].
- π = 2v/[m + v(2 − n_tot)], the exact inverse of v = m·π/[2(1−π)+π·n_tot] [1][2].
- Normal contribution fixed at 2 copies weighted (1−π) [2][4].

**Intentionally simplified:**

- Aggregation uses the median of per-variant point estimates; **consequence:** no model-fit confidence interval or joint ploidy estimate (unlike ABSOLUTE/FACETS, or `FitPurityPloidy`, ONCO-ASCAT-001) is produced.

**Implemented from CNAqc (FIN-B24 F33/F34, F62, F63):** the peak-based purity QC for simple clonal karyotypes (`analyze_peaks_common`, §4.4), the complex-karyotype (`analyze_peaks_general`) and subclonal (`analyze_peaks_subclonal`) peak analyses (§4.5), BMix mixture peaks and `n_bootstrap` (§4.6), all R-locked on seeded datasets (R random stream reproduced).

**Not implemented:**

- CCF clustering of subclonal segments (`cluster_subclonal_CCF = TRUE`, mclust `Mclust(modelNames = "E")`) — not the `analyze_peaks` default.
- Joint purity+ploidy+absolute-CN model fitting inside this VAF estimator; **users should rely on:** `OncologyAnalyzer.FitPurityPloidy` (ONCO-ASCAT-001, ASCAT runASCAT port over segment logR/BAF, after `SegmentAlleleSpecificAspcf`) for the genome-wide joint fit; ABSOLUTE [3] / FACETS [4] remain external alternatives.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | VAF-only estimator fixes m=1, n_tot=2 | Assumption | Wrong on amplified/LOH loci | accepted | ASM-02; use `EstimatePurity` for other states |
| 2 | Median aggregation | Assumption | Robust central estimate, not a fitted value | accepted | does not change the per-variant formula |
| 3 | Multiplicity bounded by n_tot in `EstimatePurity` | Validation (fixed) | m > n_tot (physically impossible; CNAqc `expectations_generalised` enumerates m ∈ 1..Major ≤ n_tot) throws `ArgumentOutOfRangeException`, e.g. (v 0.5, m 3, n_tot 2) — previously returned 1/3 | fixed | review-2026-09 B24 F7 → FIN-B24 F24 |
| 5 | `AnalyzePurityPeaks` KDE lattice | Version choice | R ≥ 4.4 `density` kernel lattice by default (densities ≈ 0.1 % different from R ≤ 4.3; `LegacyDensityCoordinates = true` reproduces R 4.3.3); convolution evaluated directly instead of FFT (≤ 1e−12 relative) | accepted | FIN-B24 F33 |
| 6 | `AnalyzePurityPeaks` degenerate inputs | Error handling | an analysed karyotype with no KDE peak in [0, 1] (R: error in `simple_peak_detector`) or no non-discarded peak throws `InvalidOperationException`; a mixture peak snapped outside histogram bins 1..100 gets `CountsPerBin = null` (R: NA / error) | accepted | FIN-B24 F34 |
| 7 | `AnalyzeSubclonalPurityPeaks` non-terminating model search | Error handling | CNAqc `evolve()` loops forever when no move stays under the ploidy cap (e.g. start 2:2 → 1:0, i.e. a 1:0 segment with `starting_state` 2:2); the port throws `InvalidOperationException` | accepted | FIN-B24 F62 |
| 8 | `AnalyzeComplexKaryotypePeaks` gate without complex karyotypes | Error handling | the `analyze_peaks` gate (a karyotype outside `karyotypes` with n > `n_min`) can open with no complex karyotype to analyse (`karyotypes` omits a simple one); CNAqc then fails (`1:nrow(NULL)`), the port throws `InvalidOperationException` | accepted | FIN-B24 F62 |
| 9 | BMix `%*%` dot products | Numerical | R evaluates `z %*% NV` with the system BLAS (OpenBLAS dgemv); the port uses R's `matprod = "internal"` long-double sum — bit-identical to R under that option, ≤ 1e−15 relative to OpenBLAS (means only; NLL/ICL and every selection identical on the 7 reference fits) | accepted | FIN-B24 F63 |
| 10 | `FitMixturePeaks` default | API default | CNAqc always fits BMix; the port keeps `FitMixturePeaks = false` by default so existing callers (VAF + karyotype only, no read counts) are unchanged — set it (with `AlternateReads`/`Depth`) to run CNAqc's default path | accepted | FIN-B24 F63 |
| 4 | Boundary π = 1 rounding in `EstimatePurity` | Numerical tolerance (fixed) | the exact CNAqc clonal peak at π = 1, v = m/n_tot, can evaluate to 1 + k·ulp (≤ 0.47·(n_tot+4)·ε for n_tot ≤ 2000); a computed π ≤ 1 + (n_tot+4)·ε (ε = 2⁻⁵²) is accepted and clamped to 1.0, e.g. (v 0.2, m 1, n_tot 5) → 1.0, (v 0.4, m 2, n_tot 5) → 1.0; larger excess still throws | fixed | review-2026-09 B24 F8 → FIN-B24 F25 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| VAF = 0 | purity = 0 | ρ = 2·0 = 0 (INV-04) |
| VAF = 0.5 (diploid) | purity = 1.0 | ρ = 2·0.5 = 1 |
| VAF > 0.5 (diploid model) | ArgumentOutOfRangeException | ρ > 1 impossible (INV-01) |
| Empty collection | ArgumentException | purity undefined [1] |
| null collection | ArgumentNullException | guard |
| n_tot < 1 or m < 1 (allele-specific) | ArgumentOutOfRangeException | formula domain |
| m > n_tot (allele-specific) | ArgumentOutOfRangeException | CNAqc m ∈ 1..Major ≤ n_tot |
| v = m/n_tot (π = 1 clonal peak) | purity = 1.0 exactly | rounding excess ≤ (n_tot+4)·ε clamped |

### 6.2 Limitations

Purity below ~0.1 approaches sequencing noise and is reported as a small value, not validated against a detection-limit model. The VAF-only estimator assumes clonality (c=1) and copy-neutral diploid heterozygosity; subclonal or amplified-segment variants must use the allele-specific overload with the correct (m, n_tot). No confidence interval, ploidy, or whole-genome-doubling handling is provided here (see ONCO-PLOIDY-001).

**LIMITATIONS:** the purity entry points answer three different questions. (1) `EstimatePurityFromVaf` / `EstimatePurityFromVAF` / `EstimatePurity` are point estimators: each variant's VAF is inverted in closed form (CNAqc `expected_vaf_fun`) and the per-variant purities are combined by their median — unweighted by read depth, with no binomial likelihood (e.g. PurBayes) and no ploidy estimation; the VAF-only overloads assume clonal heterozygous copy-neutral diploid SNVs. (2) `AnalyzePurityPeaks` / `AnalyzeComplexKaryotypePeaks` / `AnalyzeSubclonalPurityPeaks` (CNAqc `analyze_peaks`) do not estimate purity: they QC a caller-supplied purity against the VAF peaks it implies — λ = Σ weight·offset (CNAqc's "Purity correction") with a PASS/FAIL verdict for simple karyotypes, matched-peak proportions for complex karyotypes and evolution-model rankings for subclonal segments. (3) The joint purity + ploidy + absolute copy-number fit is `FitPurityPloidy` (ASCAT `runASCAT`, ONCO-ASCAT-001), with sub-clonal copy number via `FitSubclonalCopyNumber` / `CallBattenbergSubclones`. BMix mixture peaks (`FitMixturePeaks`, needs read counts) and `n_bootstrap` (`BootstrapCount`) follow CNAqc and reproduce R's random stream for a given `Seed`. Not ported from CNAqc: the non-default CCF clustering of subclonal segments (mclust).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Clonal heterozygous somatic SNV at a copy-neutral diploid locus, VAF 0.30.
double purity = OncologyAnalyzer.EstimatePurityFromVaf(0.30); // 0.60  (ρ = 2·VAF)

// Allele-specific: a 2:1 segment (n_tot=3), m=2, VAF 2/3 ⇒ fully pure tumour.
double rho = OncologyAnalyzer.EstimatePurity(new[]
{
    new OncologyAnalyzer.PurityVariant(2.0 / 3.0, Multiplicity: 2, TumorTotalCopyNumber: 3)
}); // 1.0
```

**Numerical walk-through:** CNAqc 2:1 example at π=1 [1]: m=1 ⇒ v = 1·1/[0 + 1·3] = 1/3; m=2 ⇒ v = 2/3. Inverting m=2, v=2/3, n_tot=3: π = 2·(2/3)/[2 + (2/3)(2−3)] = (4/3)/(4/3) = 1.0.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_EstimatePurity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_EstimatePurity_Tests.cs) — covers INV-01..INV-04; [OncologyAnalyzer_AnalyzePurityPeaks_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AnalyzePurityPeaks_Tests.cs) — CNAqc R-locked peak QC (16 scenarios); [OncologyAnalyzer_CnaqcGeneralSubclonalPeaks_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_CnaqcGeneralSubclonalPeaks_Tests.cs) — `analyze_peaks_general` (10 R runs) / `analyze_peaks_subclonal` (8 R runs) line-identical; [OncologyAnalyzer_CnaqcBMixBootstrap_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_CnaqcBMixBootstrap_Tests.cs) — BMix grids (7 fits, bit-identical to R `matprod = "internal"`), BMix + `n_bootstrap` pipelines (9 + 3 + 2 R runs); [StatisticsHelper_RSumDbinom_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/StatisticsHelper_RSumDbinom_Tests.cs) — R long-double `sum`, `dbinom`; [StatisticsHelper_GaussianKde_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/StatisticsHelper_GaussianKde_Tests.cs) — R `density`/`bw.nrd0`/`peakpick`
- Evidence: [ONCO-PURITY-001-Evidence.md](../../../docs/Evidence/ONCO-PURITY-001-Evidence.md)
- Related algorithms: [Variant_Allele_Frequency](../Oncology/Variant_Allele_Frequency.md)

## 8. References

1. Antonello A, Bergamin R, Calonaci N, Househam J, Milite S, Williams MJ, Anselmi F, d'Onofrio A, Sundaram V, Sosinsky A, Cross WCH, Caravagna G. 2024. Computational validation of clonal and subclonal copy number alterations from bulk tumor sequencing using CNAqc. Genome Biology 25(1):38. https://doi.org/10.1186/s13059-024-03170-5
2. CNAqc package vignette (Caravagna lab). Quality control of allele-specific copy numbers, mutations and tumour purity. https://caravagnalab.github.io/CNAqc/articles/CNAqc.html
3. Carter SL, Cibulskis K, Helman E, McKenna A, Shen H, Zack T, Laird PW, Onofrio RC, Winckler W, Weir BA, Beroukhim R, Pellman D, Levine DA, Lander ES, Meyerson M, Getz G. 2012. Absolute quantification of somatic DNA alterations in human cancer. Nature Biotechnology 30(5):413–421. https://doi.org/10.1038/nbt.2203
4. Shen R, Seshan VE. 2016. FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. Nucleic Acids Research 44(16):e131. https://doi.org/10.1093/nar/gkw520
