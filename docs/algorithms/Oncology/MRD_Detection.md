# Minimal/Molecular Residual Disease (MRD) Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-MRD-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-15 |

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

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| tumorMarkers | IEnumerable\<TumorMarker\> | required | Patient-specific tracked markers with plasma reads | non-null, non-empty |
| positivityThreshold | int | 2 | Min detected variants to call positive (τ) | ≥ 1 |
| minSupportingReads | int | 1 | Min alt reads for a marker to count as detected (r_min) | ≥ 1 |
| genomeEquivalents | int | 0 | Sequenced haploid genome equivalents n for panel Poisson p | ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Status | MrdStatus | Positive when `D ≥ τ`, else Negative |
| DetectedVariantCount | int | `D` — markers with alt reads ≥ r_min |
| TrackedVariantCount | int | `m` — panel size |
| IntegratedMutantAlleleFraction | double | IMAF = Σalt / Σtotal across loci |
| DetectionProbability | double | Panel Poisson `p = 1 − e^(−n·IMAF·m)` |

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

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Positivity threshold `τ = 2` (default), per Reinert 2019 / Signatera [1][4].
- Panel size guidance: ≤ 8 tracked SNVs compromises sensitivity at ≤ 0.1% VAF; the assay targets 16 [2].
- Poisson detection: `p = 1 − e^(−n·f·m)` (white paper Figure 2) [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DetectMRD | O(m) | O(1) | single pass over `m` markers |
| TrackVariantsOverTime | O(T·m) | O(T) | T timepoints, m markers each |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.DetectMRD(...)`: tumour-informed panel-level MRD call (≥2 detected ⇒ positive).
- `OncologyAnalyzer.TrackVariantsOverTime(...)`: longitudinal per-timepoint MRD with first-positive index.
- `OncologyAnalyzer.IsVariantDetected(...)`: per-locus presence (alt reads ≥ minSupportingReads).
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

**Intentionally simplified:**

- Per-locus detection uses a supporting-read cutoff (`r_min`, default 1); **consequence:** INVAR's exact
  per-locus trinucleotide-context background error model / GLRT is not reproduced, so the per-variant flag is
  a presence call rather than a calibrated per-locus p-value. The panel-level call and IMAF are source-exact.

**Not implemented:**

- Trinucleotide-context background error suppression and the INVAR likelihood-ratio score; **users should
  rely on:** an error-modelled per-locus caller upstream, then pass detected markers here for the panel call.
- CHIP/germline filtering of markers; **users should rely on:** ONCO-CHIP-001 (`FilterCHIP`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Per-locus detection = alt reads ≥ r_min | Assumption | Affects per-variant flag, not the ≥2 panel rule | accepted | r_min is a tunable parameter (default 1); see Evidence §Assumptions |

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
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectMRD_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMRD_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [ONCO-MRD-001-Evidence.md](../../../docs/Evidence/ONCO-MRD-001-Evidence.md)
- Related algorithms: [CtDNA_Analysis](../Oncology/CtDNA_Analysis.md)

## 8. References

1. Reinert T, Henriksen TV, Christensen E, et al. 2019. Analysis of Plasma Cell-Free DNA by Ultradeep Sequencing in Patients With Stages I to III Colorectal Cancer. *JAMA Oncology* 5(8):1124–1131. https://pubmed.ncbi.nlm.nih.gov/31070691/ (DOI:10.1001/jamaoncol.2019.0528)
2. Natera Inc. 2020. A personalized, tumor-informed approach to detect molecular residual disease with high sensitivity and specificity (Signatera analytical-validation white paper). https://www.natera.com/wp-content/uploads/2020/11/Oncology-Clinical-A-personalized-tumor-informed-approach-to-detect-molecular-residual-disease-SGN_SR_WP.pdf
3. Wan JCM, Heider K, Gale D, et al. 2020. ctDNA monitoring using patient-specific sequencing and integration of variant reads. *Science Translational Medicine* 12(548):eaaz8084. https://www.science.org/doi/10.1126/scitranslmed.aaz8084 (DOI:10.1126/scitranslmed.aaz8084)
4. Tumor-informed ctDNA MRD review (quotes the Reinert/Signatera 16-SNV, ≥2-positive rule, Table 1). PMC9265001. https://pmc.ncbi.nlm.nih.gov/articles/PMC9265001/
5. Avanzini S, et al. 2020. A mathematical model of ctDNA shedding predicts tumor detection size. *Science Advances* 6(50):eabc4308. https://doi.org/10.1126/sciadv.abc4308
