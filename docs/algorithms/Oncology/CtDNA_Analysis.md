# ctDNA Analysis (Poisson Detection, Tumour Fraction, Mean VAF)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-CTDNA-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-15 |

## 1. Overview

Circulating tumour DNA (ctDNA) analysis quantifies tumour-derived cell-free DNA (cfDNA) in plasma. This unit provides the well-defined quantitative pieces of ctDNA analysis: (a) the probability of *detecting* a variant in plasma under a Poisson sampling model given the number of input molecules and the mutant allele fraction [1][2]; (b) the limit-of-detection decision built on that model and the validated assay range [3]; (c) estimation of the ctDNA *tumour fraction* from clonal heterozygous SNV allele frequencies (tumour fraction = 2 · VAF) [4]; (d) the mean variant-allele-fraction summarization across reporters [3]; and (e) conversion of a cfDNA input mass to the haploid genome equivalents that parameterize the detection model [5][6]. The detection model is probabilistic; the tumour-fraction and mean-VAF computations are exact.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

cfDNA in plasma is a mixture of normal-derived and tumour-derived fragments. Only a small fraction is tumour-derived (median ~0.1% pre-treatment, ranging ~0.02%–3.2% by SNV/indel reporters in NSCLC) [3]. Whether a given tumour mutation is observed in plasma is governed by how many tumour molecules carrying it are actually sampled and sequenced — a counting problem in the low-burden regime [2].

### 2.2 Core Model

**Poisson detection.** With *n* sequenced haploid genome equivalents, mutant allele fraction *d*, and *k* independent tumour reporters, the number of mutant molecules sampled follows a Poisson distribution with mean λ = n·d (per reporter); the probability of observing at least one mutant molecule is [2]:

```
single reporter:  x = 1 − e^(−n·d)
k reporters:       p = 1 − e^(−n·d·k)
```

In the low-burden regime (λ < 3) detection is Poisson-limited: small changes in input or recovery shift results across the limit of detection [2].

**Mass → genome equivalents.** One haploid human genome = 3.3 pg of DNA [5]; therefore genome equivalents = mass(ng) · 1000 / 3.3 (≈ 303 GE per ng) [6].

**Tumour fraction.** For a clonal heterozygous somatic SNV at a copy-neutral diploid locus the expected VAF is v = π/2 (the m = 1, n_tot = 2 case of the CNAqc expected-VAF relation v = m·π / [2(1−π) + π·n_tot]) [4], so the tumour-derived fraction = 2·v.

**Mean VAF.** Per reporter VAF = alt reads / total reads; the ctDNA level is summarized as the mean fraction across reporters [3].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Mutant molecules are sampled independently and rarely (Poisson regime) | At very high λ the Poisson and binomial agree, so no practical impact; at extreme over-dispersion the probability is mis-estimated |
| ASM-02 | Tumour-fraction loci are clonal, heterozygous, copy-neutral diploid | VAF no longer equals TF/2; tumour fraction estimate is biased (sub-clonal or CN-altered loci) |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Detection probability p = 1 − e^(−n·d·k) ∈ [0, 1] | exponential of a non-positive exponent ∈ (0,1]; 1 minus it ∈ [0,1) [2] |
| INV-02 | p = 0 ⇔ λ = n·d·k = 0 (n = 0 or d = 0) | 1 − e⁰ = 0 [2] |
| INV-03 | p is non-decreasing in n, d, k (strictly increasing while n>0, d>0) | λ is monotone in each factor; 1 − e^(−λ) is increasing in λ [2] |
| INV-04 | tumour fraction = 2 · mean VAF, clamped to [0, 1] | v = π/2 ⇒ TF = 2v [4]; a fraction cannot exceed 1 |
| INV-05 | genome equivalents are linear in mass: GE(x ng) = x·1000/3.3, GE(0) = 0 | constant pg/genome conversion [5] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `genomeEquivalents` | `int` | required | Sequenced haploid genome equivalents n | ≥ 0 |
| `mutantAlleleFraction` | `double` | required | Mutant allele fraction d | finite, ∈ [0, 1] |
| `reporterCount` | `int` | 1 | Independent tumour reporters k | ≥ 1 |
| `minDetectionProbability` | `double` | 0.95 | Threshold p to call detected (`IsCtDnaDetected`) | ∈ (0, 1] |
| `variants` | `IEnumerable<VariantObservation>` | required | Plasma reporters (tumour read counts used) | non-null, non-empty |
| `cfDnaNanograms` | `double` | required | cfDNA input mass | finite, ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CtDnaDetectionProbability` | `double` | p = 1 − e^(−n·d·k) ∈ [0, 1] |
| `ExpectedMutantMolecules` | `double` | λ = n·d·k |
| `IsCtDnaDetected` | `bool` | λ ≥ 1 AND p ≥ threshold |
| `CalculateTumorFraction` | `double` | 2 · mean clonal het VAF ∈ [0, 1] |
| `CalculateMeanVaf` | `double` | mean of alt/total across reporters ∈ [0, 1] |
| `HaploidGenomeEquivalents` | `double` | ng · 1000 / 3.3 |

### 3.3 Preconditions and Validation

`genomeEquivalents` < 0, `mutantAlleleFraction` ∉ [0,1] or NaN, `reporterCount` < 1, and `minDetectionProbability` ∉ (0,1] throw `ArgumentOutOfRangeException`. Null variant collections throw `ArgumentNullException`; empty collections throw `ArgumentException` (the statistic is undefined). A per-variant VAF > 0.5 in `CalculateTumorFraction` throws `ArgumentOutOfRangeException` (impossible for a diploid heterozygous SNV). `cfDnaNanograms` negative or non-finite throws `ArgumentOutOfRangeException`. Read-count validation (alt ≤ total, non-negative) is delegated to the shared VAF helper.

## 4. Algorithm

### 4.1 High-Level Steps

1. Compute λ = n·d·k.
2. Detection probability: p = 1 − e^(−λ).
3. Detection decision: detected ⇔ λ ≥ 1 AND p ≥ threshold.
4. Tumour fraction: mean of per-variant VAFs (each ≤ 0.5), multiply by 2, clamp to [0,1].
5. Mean VAF: arithmetic mean of per-variant VAFs.
6. Genome equivalents: ng · 1000 / 3.3.

### 4.2 Decision Rules, Scoring, Reference Tables

| Constant | Value | Source |
|----------|-------|--------|
| pg per haploid genome | 3.3 | Devonshire et al. 2014 [5] |
| genome equivalents per ng | ≈ 303 (1000/3.3) | Alcaide et al. 2020 [6] |
| tumour fraction factor | 2 (TF = 2·VAF) | Antonello et al. 2024 [4] |
| min expected mutant molecules λ | 1 | Poisson detection model [2] |
| default detection probability | 0.95 | 95% sensitivity convention; Newman 2014 [3] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Detection probability / λ / GE | O(1) | O(1) | scalar arithmetic |
| Tumour fraction / mean VAF | O(n) | O(1) | single pass over n reporters; matches Registry O(n) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CtDnaDetectionProbability(int, double, int)`: p = 1 − e^(−n·d·k).
- `OncologyAnalyzer.ExpectedMutantMolecules(int, double, int)`: λ = n·d·k.
- `OncologyAnalyzer.IsCtDnaDetected(int, double, int, double)`: λ ≥ 1 AND p ≥ threshold.
- `OncologyAnalyzer.CalculateTumorFraction(IEnumerable<VariantObservation>)`: 2 · mean clonal het VAF, clamped.
- `OncologyAnalyzer.CalculateMeanVaf(IEnumerable<VariantObservation>)`: mean of alt/total.
- `OncologyAnalyzer.HaploidGenomeEquivalents(double)`: ng → GE.

### 5.2 Current Behavior

`CtDnaDetectionProbability` computes `1.0 - Math.Exp(-lambda)` (.NET lacks `Math.Expm1`); for the tested λ ≥ 0.01 regime this is well-conditioned. `CalculateTumorFraction` and `CalculateMeanVaf` reuse the existing private `CalculateVaf` helper (so read-count validation is shared with the somatic-calling methods). Tumour fraction is clamped to 1.0 because a fraction cannot exceed unity; per-variant VAF > 0.5 is rejected before averaging. This unit involves no substring search/matching, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Poisson detection probability p = 1 − e^(−n·d·k) and mean λ = n·d·k [2].
- Mass→genome-equivalents conversion at 3.3 pg/haploid genome [5][6].
- Tumour fraction = 2 · VAF for clonal heterozygous copy-neutral diploid SNVs [4].
- Mean variant allele fraction across reporters [3].

**Intentionally simplified:**

- Detection decision: the literature gives the probability but no single universal "called detected" threshold; the boolean `IsCtDnaDetected` adds a caller-supplied probability threshold (default 0.95) and a λ ≥ 1 physical floor. **Consequence:** the boolean (not the probability) depends on the 0.95 default; callers needing a different operating point pass their own threshold.

**Not implemented:**

- `AnalyzeFragmentSizeDistribution(bamFile)` (fragmentomics) — requires BAM-parsing infrastructure absent from this library and is a separate analysis from the quantitative ctDNA pieces in scope. **Users should rely on:** no current in-repo alternative; out of scope for ONCO-CTDNA-001.
- CHIP-background filtering and matched-tumour VCF cross-referencing — covered by separate units (ONCO-CHIP-001). **Users should rely on:** those units.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| n = 0 or d = 0 | p = 0; not detected | λ = 0 ⇒ 1 − e⁰ = 0 [2] |
| Sub-molecule input (λ < 1) | not detected even if p > 0 | cannot observe < 1 expected mutant molecule [2] |
| VAF below assay LoD (~0.025%) | low p; typically not detected | validated detection range starts at 0.025% [3] |
| Mean VAF > 0.5 in tumour-fraction | per-variant VAF > 0.5 throws | impossible for diploid het SNV [4] |
| Empty / null variant set | `ArgumentException` / `ArgumentNullException` | statistic undefined |

### 6.2 Limitations

The detection model assumes independent Poisson sampling and does not model PCR/library duplication, error-suppression (UMIs), or strand-specific recovery. Tumour-fraction estimation is valid only for clonal heterozygous copy-neutral diploid loci; sub-clonal variants or copy-number-altered loci bias the estimate. No fragmentomics (fragment-size distribution) is provided.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// n = 15,000 genome equivalents (~50 ng cfDNA), d = 0.001 (0.1% VAF), 1 reporter ⇒ λ = 15.
double p = OncologyAnalyzer.CtDnaDetectionProbability(15000, 0.001, 1); // 0.99999969...
bool detected = OncologyAnalyzer.IsCtDnaDetected(15000, 0.001);          // true

// Tumour fraction from two clonal het SNVs at VAF 0.10 and 0.20:
var variants = new[]
{
    new OncologyAnalyzer.VariantObservation("1", 100, "A", "T", 10, 100, 0, 100),
    new OncologyAnalyzer.VariantObservation("1", 200, "C", "G", 20, 100, 0, 100),
};
double tf = OncologyAnalyzer.CalculateTumorFraction(variants); // 2 * 0.15 = 0.30
```

**Numerical walk-through:** n·d = 15,000 × 0.001 = 15 expected mutant molecules (cf. Pessoa 2023) [7]; p = 1 − e^(−15) = 0.99999969409…

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CtDnaAnalysis_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CtDnaAnalysis_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [ONCO-CTDNA-001-Evidence.md](../../../docs/Evidence/ONCO-CTDNA-001-Evidence.md)
- TestSpec: [ONCO-CTDNA-001.md](../../../tests/TestSpecs/ONCO-CTDNA-001.md)
- Related algorithms: [Clonal_Subclonal_Classification](Clonal_Subclonal_Classification.md)

## 8. References

1. Newman A.M., Bratman S.V., To J., et al. 2014. An ultrasensitive method for quantitating circulating tumor DNA with broad patient coverage (CAPP-Seq). Nature Medicine 20(5):548–554. https://doi.org/10.1038/nm.3519
2. Avanzini S., Kurtz D.M., Chabon J.J., et al. 2020. A mathematical model of ctDNA shedding predicts tumor detection size. Science Advances 6(50):eabc4308. https://doi.org/10.1126/sciadv.abc4308 (detection equations corroborated verbatim via US Patent 11,085,084 B2, https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11085084)
3. Newman et al. 2014, CAPP-Seq full text. https://pmc.ncbi.nlm.nih.gov/articles/PMC4016134/
4. Antonello A., et al. 2024. CNAqc: quality control for cancer copy-number and purity. Genome Biology 25:38. https://doi.org/10.1186/s13059-024-03170-5
5. Devonshire A.S., Whale A.S., Gutteridge A., et al. 2014. Towards standardisation of cell-free DNA measurement in plasma. Anal Bioanal Chem. https://pmc.ncbi.nlm.nih.gov/articles/PMC4182654/
6. Alcaide M., Cheung M., Hillman J., et al. 2020. Evaluating the quantity, quality and size distribution of cell-free DNA by multiplex droplet digital PCR. Scientific Reports 10:12564. https://doi.org/10.1038/s41598-020-69432-x
7. Pessoa L.S., et al. 2023. Genomic approaches to cancer and minimal residual disease detection using circulating tumor DNA. https://pmc.ncbi.nlm.nih.gov/articles/PMC10314661/
