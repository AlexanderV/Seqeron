# Somatic Mutation Calling

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-SOMATIC-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Somatic mutation calling classifies each tumor variant as **somatic** (acquired by the tumor) or **germline** (inherited) by comparing the variant's allele frequency in the tumor against its allele frequency in a matched normal sample. A variant is somatic when it is present in the tumor and absent in the matched normal [1][3]. This implementation is a deterministic, rule-based classifier on observed allele fractions; it realizes the somatic-state criterion of established callers but does not reproduce their full Bayesian probability / log-odds models, which are intentionally out of scope.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A tumor genome carries both inherited germline variants and acquired somatic mutations. To isolate the somatic mutations that drive or mark the cancer, a matched normal sample (e.g. blood) from the same patient is sequenced alongside the tumor. Variants seen in both samples are germline; variants seen only in the tumor are somatic [1][3].

### 2.2 Core Model

For a variant at one site, the variant allele frequency (VAF) in a sample is `f = altReads / totalReads` (continuous allele frequency) [2]. Strelka defines the somatic state as **S = {(f_t, f_n): f_t ≠ f_n}** — the tumor allele frequency `f_t` differs from the normal allele frequency `f_n` — and restricts reported calls to a **homozygous-reference normal genotype** P(S, G_n = ref/ref | D) [1]. Mutect2 states that with no matched normal the normal likelihood ℓ_n = 1, and given a matched normal it skips variants clearly present in the germline [3].

The rule-based realization: a variant is **Somatic** when `f_t ≥ τ_t` (present in tumor) and `f_n ≤ τ_n` (absent / ref-ref in normal); **Germline** when `f_t ≥ τ_t` and `f_n > τ_n`; **NotDetected** when `f_t < τ_t`. The tumor detection threshold τ_t = 0.05 follows the whole-exome VAF limit of detection (calls ≤ 5% VAF are frequently sequencing errors) [4]; the normal absence ceiling τ_n = 0.01 is the configurable sub-percent noise band consistent with the ref/ref restriction [1][4].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Normal at homozygous reference is captured by a fixed VAF ceiling τ_n (default 1%) | If τ_n too high, true somatic calls with low normal contamination are kept; if too low, CHIP/contamination escapes as somatic |
| ASM-02 | Somatic confidence is the monotone allele-frequency separation max(0, f_t − f_n) | Not a calibrated probability; cannot be compared to caller LOD/QSS scores |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Somatic ⇔ f_t ≥ τ_t ∧ f_n ≤ τ_n | direct from the decision rule [1][4] |
| INV-02 | f_t < τ_t ⇒ NotDetected | LoD gate applied first [4] |
| INV-03 | 0 ≤ SomaticScore ≤ 1, = max(0, f_t − f_n) | both VAFs ∈ [0,1] (ASM-02) |
| INV-04 | FilterGermlineVariants returns exactly the Somatic subset | filter is defined as Status == Somatic [3] |
| INV-05 | totalReads = 0 ⇒ VAF = 0 | uncovered site treated as allele-absent [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| variants | `IEnumerable<VariantObservation>` | required | tumor variants with matched-normal read evidence | non-null |
| tumorVafThreshold | `double` | 0.05 | τ_t, min tumor VAF for presence | [0, 1] |
| normalVafThreshold | `double` | 0.01 | τ_n, max normal VAF for absence | [0, 1] |
| VariantObservation.TumorAltReads / TumorTotalReads | `int` | required | tumor alt and total read counts | 0 ≤ alt ≤ total |
| VariantObservation.NormalAltReads / NormalTotalReads | `int` | required | normal alt and total read counts | 0 ≤ alt ≤ total |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| SomaticCall.TumorVaf | `double` | f_t = tumorAlt / tumorTotal |
| SomaticCall.NormalVaf | `double` | f_n = normalAlt / normalTotal (0 if uncovered) |
| SomaticCall.Status | `SomaticStatus` | Somatic / Germline / NotDetected |
| SomaticCall.SomaticScore | `double` | max(0, f_t − f_n) when Somatic, else 0 |

### 3.3 Preconditions and Validation

`variants` must be non-null (`ArgumentNullException`). Thresholds must be in [0, 1] (`ArgumentOutOfRangeException`). Read counts must be non-negative and `alt ≤ total` (`ArgumentOutOfRangeException`). Positions are 1-based. `totalReads = 0` yields VAF 0 (tumor-only / uncovered normal). Output preserves input order and count.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs and thresholds.
2. For each variant compute `f_t = tumorAlt/tumorTotal` and `f_n = normalAlt/normalTotal` (0 if total = 0).
3. If `f_t < τ_t` → NotDetected.
4. Else if `f_n ≤ τ_n` → Somatic; else → Germline.
5. Somatic score = `max(0, f_t − f_n)` for somatic calls (0 otherwise).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Condition | Status |
|-----------|--------|
| f_t < τ_t | NotDetected |
| f_t ≥ τ_t ∧ f_n ≤ τ_n | Somatic |
| f_t ≥ τ_t ∧ f_n > τ_n | Germline |

τ_t = 0.05 [4]; τ_n = 0.01 [1][4]; diploid heterozygous fraction reference 1/2 noted by FilterMutectCalls [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CallSomaticMutations | O(n) | O(n) | n = variants; one pass, output list |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CallSomaticMutations(variants, τ_t, τ_n)`: classifies every variant.
- `OncologyAnalyzer.Classify(variant, τ_t, τ_n)`: single-variant classification.
- `OncologyAnalyzer.FilterGermlineVariants(variants, τ_t, τ_n)`: returns the somatic subset.
- `OncologyAnalyzer.CalculateSomaticScore(variant)`: separation score in [0, 1].

### 5.2 Current Behavior

Pure in-memory classification on `VariantObservation` records (no VCF parsing in this unit). This is not a search/matching operation, so the repository suffix tree is not applicable (no substring search, pattern matching, or occurrence enumeration is involved).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Somatic ⇔ present in tumor and absent in matched normal, with normal at ref/ref [1][3].
- VAF = altReads/totalReads continuous allele frequency [2].
- Tumor 5% VAF limit of detection gate [4].
- Tumor-only mode handled (no normal coverage ⇒ f_n = 0, ℓ_n = 1 analogue) [3].

**Intentionally simplified:**

- Normal absence: a fixed VAF ceiling τ_n replaces Strelka's Bayesian ref/ref genotype posterior; **consequence:** a single cutoff rather than a probability — borderline normal contamination near τ_n flips deterministically instead of by posterior.
- Somatic confidence: max(0, f_t − f_n) instead of a Bayesian somatic LOD/QSS; **consequence:** the score is a monotone surrogate, not a calibrated quality.

**Not implemented:**

- Bayesian somatic likelihood / LOD models (Strelka somaticLOD, Mutect2 TLOD/NLOD), panel-of-normals and germline-resource priors, LOH/copy-number correction; **users should rely on:** GATK Mutect2 / Strelka2 for production probabilistic calling.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | τ_n fixed normal ceiling | Assumption | borderline contamination classification | accepted | ASM-01; configurable parameter |
| 2 | Separation score | Assumption | not a calibrated probability | accepted | ASM-02 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Tumor-only (normal total = 0) | f_n = 0 ⇒ Somatic if f_t ≥ τ_t | ℓ_n = 1 with no matched normal [3] |
| Sub-5% tumor VAF | NotDetected | WES LoD [4] |
| Normal VAF just above τ_n (CHIP) | Germline | present in normal [3] |
| f_n ≥ f_t | score = 0 | no somatic separation (INV-03) |
| Empty input | empty result | trivial |

### 6.2 Limitations

Single-site allele-fraction rule only: no haplotype modeling, no indel error model, no copy-number / LOH correction, no panel-of-normals, no read-level artifact filtering. Not a substitute for a probabilistic somatic caller; thresholds are deterministic cutoffs, not calibrated posteriors.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var variants = new[]
{
    new OncologyAnalyzer.VariantObservation("chr1", 100, "A", "T", 25, 100, 0, 100),  // Somatic
    new OncologyAnalyzer.VariantObservation("chr1", 200, "C", "G", 48, 100, 50, 100), // Germline
};
var calls = OncologyAnalyzer.CallSomaticMutations(variants);
// calls[0].Status == Somatic, calls[0].SomaticScore == 0.25
// calls[1].Status == Germline, calls[1].SomaticScore == 0.0
```

**Numerical walk-through:** Variant A: f_t = 25/100 = 0.25 ≥ 0.05 and f_n = 0/100 = 0.00 ≤ 0.01 → Somatic; score = 0.25 − 0.00 = 0.25. Variant B: f_t = 0.48 ≥ 0.05 but f_n = 0.50 > 0.01 → Germline; score 0.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CallSomaticMutations_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_CallSomaticMutations_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [ONCO-SOMATIC-001-Evidence.md](../../../docs/Evidence/ONCO-SOMATIC-001-Evidence.md)

## 8. References

1. Saunders CT, Wong WSW, Swamy S, Becq J, Murray LJ, Cheetham RK. 2012. Strelka: accurate somatic small-variant calling from sequenced tumor-normal sample pairs. Bioinformatics 28(14):1811–1817. https://doi.org/10.1093/bioinformatics/bts271
2. Kim S, Scheffler K, Halpern AL, et al. 2018. Strelka2: fast and accurate calling of germline and somatic variants. Nature Methods 15:591–594. https://doi.org/10.1038/s41592-018-0051-x
3. Benjamin D, et al. / Broad Institute. 2019. GATK Mutect2 model documentation (mutect.tex). https://raw.githubusercontent.com/broadinstitute/gatk/master/docs/mutect/mutect.tex
4. Yan YH, Chen SX, Cheng LY, et al. 2021. Confirming putative variants at ≤ 5% allele frequency using allele enrichment and Sanger sequencing. Scientific Reports 11:11640. https://doi.org/10.1038/s41598-021-91142-1
