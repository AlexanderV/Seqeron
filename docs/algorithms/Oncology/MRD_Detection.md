# Minimal/Molecular Residual Disease (MRD) Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-MRD-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-23 |

## 1. Overview

Tumour-informed minimal/molecular residual disease (MRD) detection asks whether circulating tumour DNA
(ctDNA) is present in a post-treatment plasma sample, using a panel of patient-specific somatic variants
first identified in the patient's own tumour. A bespoke set of (up to 16) clonal somatic SNVs is tracked in
plasma; the sample is called **MRD-positive when at least two of the tracked variants are detected** [1][2].
This is a probabilistic detection decision, not an exact measurement: ctDNA shedding is Poisson-limited at
the molecule level [2][5], so the algorithm also reports the integrated mutant allele fraction across loci
[3] and the panel-level Poisson detection probability. It is used for post-surgical/treatment surveillance,
where MRD-positivity is strongly associated with relapse [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

After curative-intent treatment, residual tumour cells may shed ctDNA into plasma at very low allele
fractions (0.01%–0.1% VAF) [2]. A tumour-naive panel rarely covers enough of a given patient's mutations to
detect signal this low; a **tumour-informed** approach instead selects the patient's own clonal somatic SNVs
from tumour whole-exome sequencing and tracks exactly those loci, so signal can be integrated across many
informative reads [2][3].

### 2.2 Core Model

Let the panel be `m` tracked patient-specific markers. For each marker `i` let `a_i` be its plasma
alternate (mutant) supporting reads and `t_i` its total covering reads.

- **Per-locus detection:** marker `i` is *detected* when `a_i ≥ r_min` (minimum supporting reads) [3].
- **Panel call (positivity rule):** let `D = Σ_i 1[a_i ≥ r_min]` be the number of detected markers. The
  sample is **MRD-positive ⟺ D ≥ τ**, with `τ = 2` by default: "Plasma samples with at least two
  tumor-specific SNVs are defined as ctDNA-positive." [1][4]
- **Integrated mutant allele fraction (IMAF):** depth-weighted (read-pooled) plasma VAF across loci,
  `IMAF = (Σ_i a_i) / (Σ_i t_i)` [3].
- **Panel-level Poisson detection probability:** `p = 1 − e^(−n·f·m)` where `n` = haploid genome
  equivalents, `f` = ctDNA VAF (here IMAF), `m` = number of tracked mutations [2]. This reuses the
  ONCO-CTDNA-001 primitive `CtDnaDetectionProbability` (p = 1 − e^(−n·d·k)) with `k = m`.

#### INVAR-style background-subtracted, tumour-AF-weighted estimator

In addition to the read-pooled IMAF above, the analyzer implements the core, caller-reproducible part of the
INVAR pipeline [3][6]. Each tracked locus `i` carries plasma mutant reads `M_i`, depth `R_i`, tumour allele
fraction `AF_i`, and a **caller-supplied** per-locus/per-context background error rate `e_i` (from control
plasma). With per-sample ctDNA fraction `p`:

- **Per-read mixture model** [6]: the probability a read at locus `i` is mutant is
  `q_i = AF_i·(1−e_i)·p + (1−AF_i)·e_i·p + e_i·(1−p) = p·g_i + e_i·(1−p)`, with `g_i = AF_i·(1−e_i) + (1−AF_i)·e_i`.
  At `p = 0` (no ctDNA) `q_i = e_i`, the background rate; loci with higher `AF_i` and lower `e_i` carry more
  signal, so the fit is signal-to-noise (AF) weighted.
- **Log-likelihood** [6]: `logL(p) = Σ_i [ lchoose(R_i, M_i) + M_i·log(q_i) + (R_i − M_i)·log(1 − q_i) ] / n`.
- **ML ctDNA fraction `p̂`** by EM [6]: E-step `Z0_i = (1−g_i)p / ((1−g_i)p + (1−e_i)(1−p))`,
  `Z1_i = g_i·p / (g_i·p + e_i(1−p))`; M-step `p = Σ(M_i·Z1_i + (R_i − M_i)·Z0_i) / ΣR_i`
  (`initial_p = 0.01`, 200 iterations).
- **Generalised likelihood-ratio (GLRT) detection statistic** [3][6]: `LR = logL(p̂) − logL(0)`.
  Pure background ⇒ `p̂ ≈ 0`, `LR ≈ 0`; LR increases monotonically with ctDNA signal.
- **Background-subtracted, depth-weighted aggregate (IMAFv2)** [6]: per locus/context subtract background,
  `bs_i = max(0, M_i/R_i − e_i)`, then `IMAFv2 = Σ_i bs_i·R_i / Σ_i R_i`.

#### INVAR fragment-size weighting, outlier suppression, locus-noise filtering, background estimation

The full INVAR refinements are also reproduced from INVAR2 [6]:

- **Fragment-size-weighted GLRT** [6] (`calc_likelihood_ratio_with_RL`): per molecule (`DP = 1`) with mutant
  indicator `M ∈ {0,1}`, fragment length `L`, and size-likelihoods `P0(L)` (normal) and `P1(L)` (tumour), the
  wild-type and mutant read likelihoods are `L0 = (1−e)·P0·(1−p) + (1−g)·P1·p` and
  `L1 = e·P0·(1−p) + g·P1·p` (with `g = AF·(1−e)+(1−AF)·e`); `logL = Σ[M·log L1 + (1−M)·log L0]/n`; EM and the
  detection statistic `LR = logL(p̂)−logL(0)` mirror the no-size form. Tumour-derived cfDNA is shorter, so a
  short fragment has higher `P1` and is up-weighted ⇒ higher sensitivity.
- **Patient-specific outlier suppression** [6] (`repolish`): estimate the sample null ctDNA fraction
  `P_ESTIMATE = max(EM, weighted.mean(AF, TUMOUR_AF))` over background-consistent loci, then flag a locus as an
  outlier when its one-sided binomial tail `P(X ≥ mutReads | DP, P_ESTIMATE)` ≤ `α / (#loci)` (Bonferroni).
- **Locus-noise filtering + control-derived background** [6] (`createLociErrorRateTable`): per locus over
  control plasma, `e = BACKGROUND_AF = Σ(ALT_F+ALT_R)/Σ DP`; the locus passes the noise filter when alt signal
  recurs in fewer than `controlProportion` of controls AND `BACKGROUND_AF < maxBackgroundAF`. The both-strands
  filter requires `ALT_F > 0 & ALT_R > 0` (or no alt reads).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A tracked marker contributes ctDNA signal only via mutant reads above background | Background errors counted as signal inflate the detected count and IMAF |
| ASM-02 | Tracked markers are clonal, patient-specific, and independent | Subclonal/shared (e.g. CHIP) markers can bias the call (handled upstream) |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | MRD-positive ⟺ `D ≥ τ` (default 2) | direct from the calling rule [1][4] |
| INV-02 | `0 ≤ D ≤ m` | `D` counts a subset of the `m` markers |
| INV-03 | `IMAF ∈ [0, 1]` | `IMAF = Σa / Σt` with `0 ≤ a_i ≤ t_i` |
| INV-04 | `p = 1 − e^(−n·f·m) ∈ [0, 1]`, non-decreasing in `m` | exponential of a non-negative rate [2] |
| INV-05 | Longitudinal order preserved; `FirstPositiveIndex` = earliest positive (or −1) | sequential scan of timepoints |
| INV-06 | `IMAFv2 ≥ 0`; a locus with VAF ≤ background contributes 0 | `max(0, VAF − e)` then depth-weighted mean [6] |
| INV-07 | Pure-background sample ⇒ `p̂ ≈ 0` and `LR ≈ 0` (not detected) | at `p=0`, `q_i = e_i` is the MLE for background-only reads [6] |
| INV-08 | `LR` is monotone non-decreasing in injected ctDNA signal | larger true `p` raises `logL(p̂)` while `logL(0)` is fixed [6] |
| INV-09 | AF-weighted `LR` ≥ flat-pooled (mean-AF) `LR` on the same data | per-locus AF concentrates signal at high-SNR loci [3][6] |
| INV-10 | Short-tumour-fragment size-weighted `LR` ≥ no-size `LR` on the same molecules | short fragments have higher `P1` ⇒ up-weighted toward the tumour state [6] |
| INV-11 | Flat size profile (`P1 == P0`) ⇒ size-weighted `LR` == no-size `LR` | the size factor cancels in `L0`/`L1` [6] |
| INV-12 | A locus whose mutant-read count gives a binomial tail ≤ `α/n` is flagged as an outlier | one-sided `binom.test(..., "greater")` vs Bonferroni threshold [6] |
| INV-13 | Estimated `e = Σ(ALT_F+ALT_R)/Σ DP`; locus passes noise filter ⟺ `signalFrac < proportion` AND `e < maxBg` | control-pooled allele fraction and recurrence rule [6] |
| INV-14 | Both-strands pass ⟺ (`ALT_F>0` AND `ALT_R>0`) OR no alt reads | strand-bias filter [6] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| tumorMarkers | IEnumerable\<TumorMarker\> | required | Patient-specific tracked markers with plasma reads | non-null, non-empty |
| positivityThreshold | int | 2 | Min detected variants to call positive (τ) | ≥ 1 |
| minSupportingReads | int | 1 | Min alt reads for a marker to count as detected (r_min) | ≥ 1 |
| genomeEquivalents | int | 0 | Sequenced haploid genome equivalents n for panel Poisson p | ≥ 0 |
| loci (EstimateInvarSignal) | IEnumerable\<InvarLocus\> | required | Tracked loci with plasma reads, tumour AF, caller-supplied background rate | non-null; ≥1 informative locus (AF>0) |
| detectionThreshold | double | 0 | Min GLRT statistic to call ctDNA-positive | ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Status | MrdStatus | Positive when `D ≥ τ`, else Negative |
| DetectedVariantCount | int | `D` — markers with alt reads ≥ r_min |
| TrackedVariantCount | int | `m` — panel size |
| IntegratedMutantAlleleFraction | double | IMAF = Σalt / Σtotal across loci |
| DetectionProbability | double | Panel Poisson `p = 1 − e^(−n·IMAF·m)` |
| InvarSignalResult.IntegratedMutantAlleleFractionV2 | double | Background-subtracted, depth-weighted aggregate tumour fraction (IMAFv2) |
| InvarSignalResult.EstimatedTumorFraction | double | ML ctDNA fraction `p̂` (EM) under the AF-weighted mixture |
| InvarSignalResult.LikelihoodRatio | double | GLRT statistic `logL(p̂) − logL(0)` |
| InvarSignalResult.Detected | bool | `LR ≥ detectionThreshold` AND ≥1 mutant read present |
| InvarSignalResult.LocusCount | int | Number of informative loci (AF > 0) |

### 3.3 Preconditions and Validation

`tumorMarkers` null → `ArgumentNullException`; empty → `ArgumentException`. `positivityThreshold < 1`,
`minSupportingReads < 1`, or `genomeEquivalents < 0` → `ArgumentOutOfRangeException`. Read counts use
1-based positions for display only; negative read counts are clamped to 0 when summing IMAF.
`TrackVariantsOverTime` validates each timepoint panel by delegating to `DetectMRD` (null/empty panel
throws).

## 4. Algorithm

### 4.1 High-Level Steps

1. For each tracked marker, decide detection: `a_i ≥ r_min`.
2. Accumulate `D` (detected count), `Σa`, `Σt`, and panel size `m`.
3. Compute `IMAF = Σa / Σt` (0 if `Σt = 0`).
4. Compute panel Poisson `p = 1 − e^(−n·IMAF·m)`.
5. Call MRD-positive iff `D ≥ τ`.
6. Longitudinal: repeat per timepoint in order; record the earliest positive index.

INVAR signal estimation (`EstimateInvarSignal`):

1. Compute IMAFv2 over covered loci: `bs_i = max(0, M_i/R_i − e_i)`, `IMAFv2 = Σ bs_i·R_i / Σ R_i`.
2. Keep informative loci (`AF_i > 0`, `R_i > 0`); floor zero background to `1/R_i` (finite logs).
3. EM-estimate `p̂` (E/M steps above, 200 iterations from `p = 0.01`).
4. Compute `logL(0)` and `logL(p̂)`; `LR = logL(p̂) − logL(0)`.
5. Detected ⟺ `LR ≥ detectionThreshold` and at least one mutant read present.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Positivity threshold `τ = 2` (default), per Reinert 2019 / Signatera [1][4].
- Panel size guidance: ≤ 8 tracked SNVs compromises sensitivity at ≤ 0.1% VAF; the assay targets 16 [2].
- Poisson detection: `p = 1 − e^(−n·f·m)` (white paper Figure 2) [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DetectMRD | O(m) | O(1) | single pass over `m` markers |
| TrackVariantsOverTime | O(T·m) | O(T) | T timepoints, m markers each |
| EstimateInvarSignal | O(m·I) | O(m) | I = EM iterations (200); m informative loci |
| IntegratedMutantAlleleFractionV2 | O(m) | O(1) | single pass over loci |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.DetectMRD(...)`: tumour-informed panel-level MRD call (≥2 detected ⇒ positive).
- `OncologyAnalyzer.TrackVariantsOverTime(...)`: longitudinal per-timepoint MRD with first-positive index.
- `OncologyAnalyzer.IsVariantDetected(...)`: per-locus presence (alt reads ≥ minSupportingReads).
- `OncologyAnalyzer.EstimateInvarSignal(...)`: INVAR-style background-subtracted, AF-weighted ctDNA estimate (IMAFv2, ML `p̂`, GLRT, detection call).
- `OncologyAnalyzer.IntegratedMutantAlleleFractionV2(...)`: background-subtracted, depth-weighted aggregate (IMAFv2).
- `OncologyAnalyzer.EstimateInvarSignalWithSize(...)`: INVAR fragment-size-weighted GLRT (with-RL) over per-molecule `InvarMolecule` + `FragmentSizeProfile`.
- `OncologyAnalyzer.SuppressOutlierLoci(...)`: INVAR patient-specific outlier suppression (Bonferroni one-sided binomial test).
- `OncologyAnalyzer.EstimateLocusBackground(...)`: control-derived per-locus background error rate + locus-noise verdict.
- `OncologyAnalyzer.PassesBothStrandsFilter(...)`: INVAR both-strands filter.
- Reuses `OncologyAnalyzer.CtDnaDetectionProbability(...)` (ONCO-CTDNA-001) for the panel Poisson `p`.

### 5.2 Current Behavior

The panel-level positivity rule and IMAF are computed in a single pass. The Poisson detection probability is
delegated to the existing ctDNA primitive (no duplicated formula). No substring/pattern search is involved,
so the repository suffix tree is **not applicable** (markers are matched positionally by the caller, not by
sequence search).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Panel positivity rule `D ≥ 2` ⇒ MRD-positive [1][4].
- IMAF = depth-weighted plasma VAF across loci [3].
- Panel Poisson detection `p = 1 − e^(−n·f·m)` [2].
- INVAR per-locus background subtraction with caller-supplied `e_i`, and the background-subtracted,
  depth-weighted aggregate IMAFv2 = `weighted.mean(max(0, VAF − e), depth)` [6] (`calculateIMAFv2`).
- INVAR AF-weighted mixture model `q = p·g + e(1−p)`, EM ML estimate of the ctDNA fraction `p̂`, and the
  generalised likelihood-ratio statistic `LR = logL(p̂) − logL(0)` [3][6] (`calc_log_likelihood`,
  `estimate_p_EM`, `calc_likelihood_ratio`, no-size variant).
- INVAR fragment-size-weighted GLRT `L0 = (1−e)·P0·(1−p) + (1−g)·P1·p`, `L1 = e·P0·(1−p) + g·P1·p` with the
  with-RL EM and `LR = logL(p̂)−logL(0)` [6] (`calc_log_likelihood_with_RL`, `estimate_p_EM_with_RL`,
  `calc_likelihood_ratio_with_RL`).
- INVAR patient-specific outlier suppression: `P_THRESHOLD = α/n`,
  `P_ESTIMATE = max(EM, weighted.mean(AF, TUMOUR_AF))`, one-sided binomial tail test [6] (`repolish`).
- INVAR control-derived background error `BACKGROUND_AF = Σ(ALT_F+ALT_R)/Σ DP`, locus-noise filter
  `(N_with_signal/N) < proportion AND BACKGROUND_AF < maxBg`, and both-strands filter [6]
  (`createLociErrorRateTable`).

**Intentionally simplified:**

- The per-variant *detected* flag of `DetectMRD` uses a supporting-read cutoff (`r_min`, default 1);
  **consequence:** the legacy ≥2-of-N panel call is a presence call. The calibrated per-locus signal is now
  available via `EstimateInvarSignal` (GLRT). The legacy panel call and read-pooled IMAF are unchanged.
- The fragment-size likelihood uses the discrete empirical proportion `COUNT/TOTAL` per length bin rather than
  INVAR2's weighted KDE (`estimate_real_length_probability`, bandwidth `adjust = 0.03`); **consequence:** on
  sparsely-populated length bins the exact size weight differs slightly from the smoothed estimate, but the
  with-RL GLRT and the short-fragment up-weighting are exact.

**Not implemented:**

- CHIP/germline filtering of markers; **users should rely on:** ONCO-CHIP-001 (`FilterCHIP`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Per-locus detection = alt reads ≥ r_min | Assumption | Affects per-variant flag, not the ≥2 panel rule | accepted | r_min is a tunable parameter (default 1); see Evidence §Assumptions |
| 2 | INVAR background rate `e_i` may be caller-supplied OR estimated from controls | Assumption | GLRT quality depends on the background model | accepted | `EstimateLocusBackground` now derives `e_i` from control-plasma reads (`BACKGROUND_AF`); `InvarLocus.BackgroundErrorRate` still accepts a caller value |
| 3 | Fragment-size likelihood uses empirical `COUNT/TOTAL` per length bin, not INVAR2's weighted KDE | Assumption | Slightly different weight on sparse length bins | accepted | with-RL GLRT/EM are exact; KDE bandwidth smoothing (`adjust = 0.03`) not reproduced (Evidence §Assumptions) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Exactly 1 detected | MRD-negative | below τ = 2 [1][4] |
| 0 detected | MRD-negative | no signal |
| All total reads = 0 | IMAF = 0, p = 0 | no informative reads |
| Empty panel | ArgumentException | nothing to interrogate |
| Custom τ = 1 | single detected ⇒ positive | parameterized threshold |

### 6.2 Limitations

Detection sensitivity is not modelled per-locus from sequencing error; with very few tracked markers
(≤ 8) low-VAF MRD may be missed [2]. The algorithm assumes markers are clonal and patient-specific;
CHIP/germline contamination must be removed upstream. The Poisson `p` is informative only when
`genomeEquivalents` (n) reflects the sequenced input.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var panel = new[]
{
    new OncologyAnalyzer.TumorMarker("1", 100, "A", "T", 3, 200),
    new OncologyAnalyzer.TumorMarker("2", 200, "C", "G", 1, 150),
    new OncologyAnalyzer.TumorMarker("3", 300, "G", "A", 0, 180),
};
OncologyAnalyzer.MrdResult mrd = OncologyAnalyzer.DetectMRD(panel);
// DetectedVariantCount = 2  (loci 1 and 2 have alt reads), Status = Positive,
// IMAF = (3+1+0)/(200+150+180) = 4/530 ≈ 0.0075472.

// INVAR-style background-subtracted, AF-weighted estimate (caller supplies background e per locus):
var loci = new[]
{
    new OncologyAnalyzer.InvarLocus(5, 1000, 0.4, 0.001),  // alt, total, tumourAF, background
    new OncologyAnalyzer.InvarLocus(5, 1000, 0.4, 0.001),
    // ... many loci ...
};
OncologyAnalyzer.InvarSignalResult inv = OncologyAnalyzer.EstimateInvarSignal(loci);
// EstimatedTumorFraction p̂ ≈ injected ctDNA fraction; LikelihoodRatio grows with signal; Detected when LR ≥ threshold.
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectMRD_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMRD_Tests.cs) — covers `INV-01`..`INV-14`
- Evidence: [ONCO-MRD-001-Evidence.md](../../../docs/Evidence/ONCO-MRD-001-Evidence.md)
- Related algorithms: [CtDNA_Analysis](../Oncology/CtDNA_Analysis.md)

## 8. References

1. Reinert T, Henriksen TV, Christensen E, et al. 2019. Analysis of Plasma Cell-Free DNA by Ultradeep Sequencing in Patients With Stages I to III Colorectal Cancer. *JAMA Oncology* 5(8):1124–1131. https://pubmed.ncbi.nlm.nih.gov/31070691/ (DOI:10.1001/jamaoncol.2019.0528)
2. Natera Inc. 2020. A personalized, tumor-informed approach to detect molecular residual disease with high sensitivity and specificity (Signatera analytical-validation white paper). https://www.natera.com/wp-content/uploads/2020/11/Oncology-Clinical-A-personalized-tumor-informed-approach-to-detect-molecular-residual-disease-SGN_SR_WP.pdf
3. Wan JCM, Heider K, Gale D, et al. 2020. ctDNA monitoring using patient-specific sequencing and integration of variant reads. *Science Translational Medicine* 12(548):eaaz8084. https://www.science.org/doi/10.1126/scitranslmed.aaz8084 (DOI:10.1126/scitranslmed.aaz8084)
4. Tumor-informed ctDNA MRD review (quotes the Reinert/Signatera 16-SNV, ≥2-positive rule, Table 1). PMC9265001. https://pmc.ncbi.nlm.nih.gov/articles/PMC9265001/
5. Avanzini S, et al. 2020. A mathematical model of ctDNA shedding predicts tumor detection size. *Science Advances* 6(50):eabc4308. https://doi.org/10.1126/sciadv.abc4308
6. Rosenfeld lab (nrlab-CRUK). INVAR2 — restructured INVAR pipeline (reference implementation of [3]). `R/shared/detectionFunctions.R`, `R/4_detection/generalisedLikelihoodRatioTest.R`, `R/3_outlier_suppression/outlierSuppression.R`, `R/3_outlier_suppression/sizeCharacterisation.R`, `R/1_parse/onTargetErrorRatesAndFilter.R`. https://github.com/nrlab-CRUK/INVAR2
7. Lanczos C. 1964. A precision approximation of the gamma function. *J. SIAM Numer. Anal.* 1(1):86–96. https://doi.org/10.1137/0701008
