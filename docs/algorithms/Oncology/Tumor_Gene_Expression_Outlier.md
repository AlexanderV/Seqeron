# Tumor Gene Expression Outlier and Signature Score

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-EXPR-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-15 |

## 1. Overview

This unit provides the oncology-specific layer of tumor gene-expression analysis: given a sample's per-gene expression and caller-supplied reference cohorts, it computes the expression **z-score** of each gene and flags **over-/under-expression outliers**, and it combines per-gene z-scores into a **gene-signature (combined z-score) activity** for a caller-defined signature. It is specification-driven (exact arithmetic, not heuristic): the z-score and the combined z-score follow the cBioPortal mRNA normalization spec [1][2] and the Lee et al. (2008) pathway-activity definition [4]. TPM/FPKM quantification and differential expression already live in `TranscriptomeAnalyzer`; this layer assumes expression values are already on a normalization scale on which a z-score is meaningful.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Cancer genes are frequently dysregulated only in a subset of tumors. A standard way to detect such single-sample over-/under-expression is to standardize a gene's expression against a reference distribution and flag large deviations [1][3]. Multi-gene signatures (e.g. prognostic or pathway panels) summarize coordinated expression of a defined gene set into one activity score [4].

### 2.2 Core Model

**Per-gene z-score** [1][2]:

```
z = (r − μ) / σ
```

where `r` is the sample expression value, and `μ`, `σ` are the mean and standard deviation of the gene's reference cohort. The reference implementation `NormalizeExpressionLevels.java` computes `σ` as the **sample** standard deviation with divisor (n − 1) [2]:

```
μ = (Σ rᵢ) / n
σ = sqrt( Σ(rᵢ − μ)² / (n − 1) )
```

**Outlier rule** [3]: a gene is an outlier iff `z > +t` (overexpressed) or `z < −t` (underexpressed); the default `t = 2`. The comparison is strict — `|z| = t` is not an outlier.

**Combined z-score (signature / pathway activity)** [4]: for a signature with k member genes whose per-gene z-scores are z₁…z_k,

```
a = (Σᵢ zᵢ) / √k
```

The √k denominator (rather than k) stabilizes the variance of the mean [4][5].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Reference cohort approximates the gene's expression distribution | z-scores miscalibrated; outlier calls biased if the cohort is unrepresentative [3] |
| ASM-02 | Inputs are on a scale where a z-score is meaningful (e.g. log-TPM, microarray intensity) | z is computed regardless, but interpretation as "standard deviations from mean" is weakened |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | z = 0 when value = cohort mean | z = (μ − μ)/σ = 0 [1] |
| INV-02 | z is monotone increasing in the value for a fixed cohort | z is linear in `r`, σ > 0 [1] |
| INV-03 | z(2μ − r) = −z(r) (reflection about the mean negates z) | linearity of z [1] |
| INV-04 | Outlier iff z > +t (Over) or z < −t (Under), strict | cBioPortal FAQ ">2 or <-2" [3] |
| INV-05 | Signature score of k equal z-scores all = c is c·√k | a = (k·c)/√k = c√k [4] |
| INV-06 | Signature score with k = 1 equals the single z-score | a = z/√1 = z [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| value | double | required | Sample expression value `r` | on a normalization scale |
| referenceCohort | IReadOnlyList\<double\> | required | Reference values for the gene | n ≥ 2; non-zero spread |
| sampleExpression | IReadOnlyDictionary\<string,double\> | required | Per-gene sample values | non-null |
| referenceCohorts | IReadOnlyDictionary\<string,IReadOnlyList\<double\>\> | required | Per-gene reference cohorts | must contain every sampled gene |
| threshold | double | 2.0 | Absolute z-score outlier cutoff | > 0 |
| memberZScores | IReadOnlyList\<double\> | required | Per-gene z-scores of a signature | k ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (z-score) | double | z = (value − μ)/σ |
| ExpressionOutlier.Gene | string | Outlier gene id |
| ExpressionOutlier.ZScore | double | Its z-score |
| ExpressionOutlier.Direction | ExpressionDirection | Over (z > +t) or Under (z < −t) |
| (signature score) | double | a = (Σ z)/√k |

### 3.3 Preconditions and Validation

Null cohorts/dictionaries → `ArgumentNullException`. Reference cohort with n < 2 → `ArgumentException` (sample SD undefined). Reference cohort with σ = 0 → `ArgumentException` (mirrors the reference implementation's fatal error [2]). Non-positive threshold → `ArgumentOutOfRangeException`. A sampled gene with no reference cohort → `ArgumentException`. Empty signature (k = 0) → `ArgumentException`. Identifiers are case-sensitive (dictionary keys). No alphabet/sequence input is involved.

## 4. Algorithm

### 4.1 High-Level Steps

1. **z-score:** compute cohort mean μ and sample SD σ (divisor n − 1); return (value − μ)/σ.
2. **Outliers:** for each sampled gene, look up its reference cohort, compute z, and emit it as Over if z > +t or Under if z < −t.
3. **Signature score:** sum the member-gene z-scores and divide by √k.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Quantity | Value | Source |
|----------|-------|--------|
| SD divisor | n − 1 (sample SD) | NormalizeExpressionLevels.java `std()` [2] |
| Default outlier threshold | 2.0, strict | cBioPortal FAQ [3] |
| Signature denominator | √k | Lee et al. (2008) [4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateExpressionZScore | O(n) | O(1) | n = cohort size; two passes (mean, SD) |
| IdentifyOutlierGenes | O(g·n) | O(g) | g = sampled genes, n = cohort size |
| CalculateSignatureScore | O(k) | O(1) | k = signature size |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CalculateExpressionZScore(double, IReadOnlyList<double>)`: per-gene z-score.
- `OncologyAnalyzer.IdentifyOutlierGenes(IReadOnlyDictionary<string,double>, IReadOnlyDictionary<string,IReadOnlyList<double>>, double)`: outlier detection.
- `OncologyAnalyzer.CalculateSignatureScore(IReadOnlyList<double>)`: combined z-score signature activity.

### 5.2 Current Behavior

Reference cohorts and the signature gene set are caller-supplied; the unit does not bundle any cohort or signature (hence Framework status). Outliers are returned in the iteration order of the sample dictionary. No search/matching is performed, so the repository suffix tree is **not applicable** to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- z = (r − μ)/σ with σ the sample SD (divisor n − 1) [1][2].
- Strict outlier rule z > +t / z < −t, default t = 2 [3].
- Combined z-score a = (Σ z)/√k [4].
- σ = 0 / no-spread reference treated as an error, mirroring the reference implementation [2].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Diploid-only base-population selection (Method 1 of cBioPortal, which restricts the reference to copy-number-neutral samples); **users should rely on:** supplying the desired reference cohort directly (this layer is reference-agnostic) [1].
- Rank-based single-sample signature methods (ssGSEA, singscore); **users should rely on:** the combined z-score here or external tools [5].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Caller supplies reference/signature | Assumption | scope of inputs | accepted | ASM-01 |
| 2 | Input on z-meaningful scale | Assumption | interpretation | accepted | ASM-02 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| value = cohort mean | z = 0 | INV-01 |
| z = ±threshold exactly | NOT an outlier | strict rule [3] |
| zero-spread cohort | ArgumentException | reference impl aborts [2] |
| cohort size 1 | ArgumentException | sample SD undefined [2] |
| signature k = 1 | a = z₁ | INV-06 |
| empty signature | ArgumentException | a undefined [4] |

### 6.2 Limitations

The z-score is only as good as the supplied reference cohort; with tumor-only cohorts the ±2 default "oftentimes [is] not meaningful" [3]. No batch-effect correction, no matched-normal modeling, and no copy-number-aware base-population selection are performed (caller responsibilities).

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** cohort {2, 2, 4, 6, 6}: μ = 4; Σ(rᵢ−μ)² = 4+4+0+4+4 = 16; sample variance = 16/(5−1) = 4; σ = 2. For value 10: z = (10−4)/2 = **3.0** (Over). For value −1: z = −2.5 (Under). For value 8: z = 2.0 (boundary → not an outlier). Signature z = {3, 1, −1, 1}: a = (3+1−1+1)/√4 = 4/2 = **2.0** [1][4].

**API usage example:**

```csharp
double z = OncologyAnalyzer.CalculateExpressionZScore(10.0, new[] { 2.0, 2.0, 4.0, 6.0, 6.0 }); // 3.0
var outliers = OncologyAnalyzer.IdentifyOutlierGenes(sample, cohorts); // genes with |z| > 2
double activity = OncologyAnalyzer.CalculateSignatureScore(new[] { 3.0, 1.0, -1.0, 1.0 }); // 2.0
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_IdentifyOutlierGenes_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_IdentifyOutlierGenes_Tests.cs) — covers INV-01..INV-06
- Evidence: [ONCO-EXPR-001-Evidence.md](../../../docs/Evidence/ONCO-EXPR-001-Evidence.md)
- Related algorithms: [Tumor_Heterogeneity_Analysis](./Tumor_Heterogeneity_Analysis.md)

## 8. References

1. cBioPortal. 2024. mRNA Expression Z-Scores normalization specification. https://docs.cbioportal.org/z-score-normalization-script/
2. cBioPortal core. NormalizeExpressionLevels.java (reference implementation; `std()` divisor n−1; σ=0 fatal error). https://github.com/cBioPortal/cbioportal-core/blob/master/src/main/java/org/mskcc/cbio/portal/scripts/NormalizeExpressionLevels.java
3. cBioPortal. FAQ — default mRNA z-score threshold ±2 (strict >2 / <−2). https://docs.cbioportal.org/user-guide/faq/
4. Lee E, Chuang H-Y, Kim J-W, Ideker T, Lee D. 2008. Inferring Pathway Activity toward Precise Disease Classification. PLoS Comput Biol 4(11):e1000217. https://doi.org/10.1371/journal.pcbi.1000217
5. Hänzelmann S, Castelo R, Guinney J. 2013. GSVA: gene set variation analysis (combined z-score method, Lee et al. 2008). https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html
