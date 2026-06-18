# Differential Expression

| Field | Value |
|-------|-------|
| Algorithm Group | Transcriptome Analysis |
| Test Unit ID | TRANS-DIFF-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Differential expression analysis identifies genes whose expression differs between two conditions (for example treatment versus control). For each gene it reports a **log2 fold change** (the log2 ratio of mean expression between the two conditions) and a **p-value** from a two-sample test, then corrects the p-values across all genes for multiple testing and flags genes that are both biologically and statistically significant. This implementation uses the classical replicate-level estimator: mean log2 ratio for effect size, Welch's unequal-variance t-test for significance, and the Benjamini-Hochberg step-up procedure for false-discovery-rate control. It is a specification-driven statistical procedure (exact given the input replicates), not a probabilistic model fit. [1][2][3][5]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

RNA-seq and microarray studies measure expression of thousands of genes across replicated samples. A gene is "differentially expressed" between two conditions when its expression changes by a meaningful magnitude and the change is unlikely under the null hypothesis of no difference. Because thousands of hypotheses are tested simultaneously, raw p-values must be adjusted to control the false discovery rate (FDR). [1][5]

### 2.2 Core Model

**Log2 fold change.** For replicate expression vectors of condition 1 (control) and condition 2 (treatment) with means `m1` and `m2`:

```
log2FC = log2((m2 + c) / (m1 + c))
```

a positive value means higher expression in the treatment (numerator). The unregularized definition is `log2(m2/m1)`; a pseudocount `c` is added to both means to keep the ratio finite when a mean is zero. [1][5]

**Welch's two-sample t-test** (unequal variances) on the replicates, with unbiased (N-1) sample variances `s1², s2²` and sizes `N1, N2`:

```
t  = (m2 - m1) / sqrt(s1²/N1 + s2²/N2)
ν  = (s1²/N1 + s2²/N2)² / [ s1⁴/(N1²(N1-1)) + s2⁴/(N2²(N2-1)) ]
```

The two-sided p-value is the exact Student's t tail via the regularized incomplete beta function:

```
p = I_{ν/(ν+t²)}(ν/2, 1/2)
```

[3][4]

**Benjamini-Hochberg FDR.** With `m` p-values ordered ascending `p(1) ≤ … ≤ p(m)`, the BH-adjusted p-value (matching R's `p.adjust(method="BH")`) is the cumulative minimum, from the largest p-value downward, of `m·p(i)/i`, clamped to 1. The decision form: reject the largest `k` with `p(k) ≤ (k/m)α`. [2]

**Significance gate.** A gene is differentially expressed when **both** `|log2FC| ≥ threshold` **and** `adjusted p-value < α`. [1][5]

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `CalculateFoldChange(a,b) = −CalculateFoldChange(b,a)` | log2((m2+c)/(m1+c)) = −log2((m1+c)/(m2+c)) |
| INV-02 | Equal group means ⇒ log2FC = 0 and p = 1 ⇒ not significant | ratio = 1 ⇒ log2 = 0; t = 0 ⇒ p = I_x(ν/2,½) at x=1 = 1 |
| INV-03 | Adjusted p ≥ raw p and ≤ 1 | BH multiplies by m/rank ≥ 1 and clamps to 1 [2] |
| INV-04 | Adjusted p-values are monotone non-decreasing in raw-p order | BH cumulative minimum from largest p down [2] |
| INV-05 | Significant ⇔ `|log2FC| ≥ threshold` AND `adjusted p < α` | two-criterion gate [1][5] |

### 2.5 Comparison with Related Methods

| Aspect | This unit (classical t-test + BH) | DESeq2 (Love et al. 2014) |
|--------|-----------------------------------|---------------------------|
| Effect size | mean log2 ratio (with pseudocount) | shrunken GLM log2 coefficient |
| Significance | Welch t-test, exact t tail | Wald test on shrunken LFC / z-statistic |
| Variance model | per-gene sample variance | dispersion shrinkage across genes |
| Multiple testing | Benjamini-Hochberg | Benjamini-Hochberg |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `genes` | `IEnumerable<(string GeneId, IReadOnlyList<double> Condition1, IReadOnlyList<double> Condition2)>` | required | Per-gene replicate expression for the two conditions | null → empty result; per-group N ≥ 2 to be testable |
| `alpha` | `double` | 0.05 | Adjusted-p (FDR) significance threshold | 0 < α ≤ 1 |
| `log2FoldChangeThreshold` | `double` | 1.0 | Minimum |log2 fold change| for significance | ≥ 0 |
| `expression1` | `IReadOnlyList<double>` | required | Replicates of condition 1 (control / reference) | null/empty → 0 |
| `expression2` | `IReadOnlyList<double>` | required | Replicates of condition 2 (treatment, numerator) | null/empty → 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `GeneId` | string | Gene identifier |
| `Log2FoldChange` | double | log2((m2+c)/(m1+c)); positive = up in condition 2 |
| `PValue` | double | Raw two-sided Welch t-test p-value |
| `AdjustedPValue` | double | Benjamini-Hochberg FDR-adjusted p-value |
| `IsSignificant` | bool | True iff |log2FC| ≥ threshold AND adjusted p < α |
| `Regulation` | string | "Upregulated" / "Downregulated" / "Unchanged" by log2FC sign |

### 3.3 Preconditions and Validation

`null` gene enumerable → empty result. `null`/empty expression lists → fold change 0. A group with fewer than 2 replicates cannot form an unbiased variance, so its p-value is set to 1 (not testable). Zero standard error: identical means → p = 1; separated means → p = 0. Inputs are treated as already-normalized expression values (no internal library-size normalization). No exceptions are thrown for degenerate input; degenerate cases resolve to the documented conventions.

## 4. Algorithm

### 4.1 High-Level Steps

1. For each gene compute `log2FC = log2((mean2+c)/(mean1+c))`.
2. For each gene compute the raw two-sided Welch t-test p-value (exact Student's t tail).
3. Adjust all raw p-values with the Benjamini-Hochberg step-up procedure.
4. Flag a gene significant iff `|log2FC| ≥ threshold` AND `adjusted p < α`; label regulation by the sign of log2FC.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Pseudocount `c = 1` for the fold-change ratio (degenerate-input regularization). [1][5]
- Welch-Satterthwaite degrees of freedom and the regularized incomplete beta `I_x(a,b)` (Lentz continued fraction) for the exact t tail. [3][4]
- BH adjustment = cumulative minimum of `m·p(i)/i` from largest rank down, clamped to 1 (R `p.adjust` BH). [2]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindDifferentiallyExpressed` | O(g·s + g·log g) | O(g) | g genes, s replicates per gene; the g·log g term is the BH sort |
| `CalculateFoldChange` | O(s) | O(1) | two mean passes |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [TranscriptomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs)

- `TranscriptomeAnalyzer.CalculateFoldChange(expression1, expression2)`: log2 fold change between two conditions.
- `TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes, alpha, log2FoldChangeThreshold)`: per-gene fold change + Welch t-test p-value + BH adjustment + two-criterion significance gate.

### 5.2 Current Behavior

Expression inputs are taken as already-normalized values; no internal library-size or TPM normalization is applied. The two-sided p-value uses the exact Student's t tail (regularized incomplete beta), not a normal approximation. This unit does not use the repository suffix tree: there is no substring search, pattern matching, or occurrence enumeration — the computation is purely numeric over expression vectors, so the suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- log2 fold change as the log2 ratio of mean expression between treatment and control. [1][5]
- Welch's unequal-variance t-statistic with unbiased (N-1) variances and Welch-Satterthwaite df. [3]
- Exact two-sided Student's t p-value `I_{ν/(ν+t²)}(ν/2, ½)`. [4]
- Benjamini-Hochberg step-up adjusted p-values (R `p.adjust` BH algorithm). [2]
- Two-criterion DE gate: |log2FC| ≥ threshold AND adjusted p < α. [1][5]

**Intentionally simplified:**

- Effect size uses the mean log2 ratio with pseudocount instead of DESeq2's dispersion-shrunken GLM coefficient; **consequence:** fold-change magnitudes for low-count genes are not shrunk toward zero as DESeq2 would, and a pseudocount of 1 is added. [1]

**Not implemented:**

- DESeq2 dispersion estimation / LFC shrinkage and the Wald-test z-statistic; **users should rely on:** DESeq2/edgeR for production count-model DE; no in-repo alternative.
- Internal library-size / TPM normalization of raw counts; **users should rely on:** `TranscriptomeAnalyzer.CalculateTPM` / `QuantileNormalize` upstream.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Fold-change pseudocount = 1 | Assumption | Only affects genes with a near-zero mean; never changes sign or non-degenerate ordering | accepted | Evidence Assumption 1 |
| 2 | <2 replicates ⇒ p = 1 | Assumption | Single-replicate genes never flagged significant | accepted | Welch precondition; Evidence Assumption 2 |
| 3 | se = 0 ⇒ p = 1 (equal means) / 0 (unequal) | Assumption | Degenerate zero-variance groups | accepted | limit of t-statistic; Evidence Assumption 3 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty gene enumerable / null | Empty result | No genes to test |
| Identical groups | log2FC = 0, p = 1, not significant | ratio = 1, t = 0 (INV-02) |
| Group with <2 replicates | p = 1, not significant | variance undefined (Welch precondition) |
| Zero mean expression | finite log2FC via pseudocount | regularized ratio |
| Strong p but |log2FC| below threshold | not significant | two-criterion gate (INV-05) |

### 6.2 Limitations

Assumes already-normalized expression and approximately normal replicate values (t-test assumption); with very few replicates the t-test has low power. Does not model count-data dispersion, batch effects, or covariates. Pseudocount choice affects only near-zero-mean genes. For production count-based DE on raw RNA-seq counts, dedicated tools (DESeq2, edgeR) are more appropriate.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** control = {1,2,3}, treatment = {7,8,9}.
means 2 and 8; sample variances 1 and 1; se = √(1/3+1/3) = 0.8164966; t = (8−2)/0.8164966 = 7.3484692; Welch df ν = (2/3)²/(2·(1/3)²/2) = 4; two-sided p = `I_{4/(4+t²)}(2, ½)` = 0.0018262607 (cross-checked against SciPy `ttest_ind(equal_var=False)`). With three additional genes giving raw p-values (0.001, 0.4, 0.5, 0.9), BH adjustment yields (0.004, 0.66667, 0.66667, 0.9). [2][3][4]

**API usage example:**

```csharp
var genes = new (string, IReadOnlyList<double>, IReadOnlyList<double>)[]
{
    ("G1", new double[] {1,2,3}, new double[] {7,8,9}),
};
var results = TranscriptomeAnalyzer
    .FindDifferentiallyExpressed(genes, alpha: 0.05, log2FoldChangeThreshold: 1.0)
    .ToList();
// results[0].Log2FoldChange > 0  → "Upregulated" in condition 2
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [TranscriptomeAnalyzer_DifferentialExpression_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_DifferentialExpression_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [TRANS-DIFF-001-Evidence.md](../../../docs/Evidence/TRANS-DIFF-001-Evidence.md)
- Related algorithms: [Expression_Quantification](./Expression_Quantification.md)

## 8. References

1. Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15:550. https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/
2. Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57(1):289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x (procedure via https://en.wikipedia.org/wiki/False_discovery_rate ; algorithm via R p.adjust https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
3. Welch BL. 1947. The generalization of "Student's" problem when several different population variances are involved. Biometrika 34(1–2):28–35. https://en.wikipedia.org/wiki/Welch%27s_t-test
4. Student's t-distribution — cumulative distribution function via regularized incomplete beta function. https://en.wikipedia.org/wiki/Student%27s_t-distribution
5. Science Park Study Group. Introduction to RNA-seq — 06 Differential expression analysis. https://scienceparkstudygroup.github.io/rna-seq-lesson/06-differential-analysis/index.html
