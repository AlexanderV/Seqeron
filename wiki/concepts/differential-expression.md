---
type: concept
title: "Differential expression analysis (log2 fold change + Welch's t-test + Benjamini-Hochberg FDR)"
tags: [transcriptome, algorithm]
mcp_tools:
  - differential_expression
sources:
  - docs/algorithms/Transcriptome/Differential_Expression.md
  - docs/Evidence/TRANS-DIFF-001-Evidence.md
source_commit: aba976debe8dc3f316f6eb6d0cb3eb337f648314
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-diff-001-evidence
      evidence: "Test Unit ID: TRANS-DIFF-001, Algorithm: Differential Expression (log2 fold change, Welch's t-test, Benjamini-Hochberg FDR); Methods TranscriptomeAnalyzer.CalculateFoldChange / FindDifferentiallyExpressed"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:significant-taxa-detection
      source: trans-diff-001-evidence
      evidence: "Both are per-feature two-group significance tests, but DE uses a parametric Welch (unequal-variance) t-test on log2-fold-change effect sizes AND applies Benjamini-Hochberg FDR across features, whereas significant-taxa uses the non-parametric Mann-Whitney U rank-sum with no built-in multiple-testing correction (the DMR concept notes BH is the caller's responsibility there)."
      confidence: medium
      status: current
---

# Differential expression analysis (log2 fold change + Welch's t-test + Benjamini-Hochberg FDR)

Two-group **RNA-seq differential expression (DE)**: for each gene, decide whether its expression
differs between a **treatment** group and a **control** group. The unit combines three
well-established pieces — a **log2 fold-change** effect size, a **Welch (unequal-variance) two-sample
t-test** raw p-value, and a **Benjamini-Hochberg (BH) FDR** adjustment across all genes — and calls a
gene DE only when **both** an effect-size threshold and an adjusted-p threshold are met. This is the
**first ingested unit of the Transcriptome / RNA-seq family**. Validated under test unit
**TRANS-DIFF-001**; the record is [[trans-diff-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the artifact pattern.

Impl `TranscriptomeAnalyzer` (`Seqeron.Genomics.Annotation`, the Annotation server):
`CalculateFoldChange(...)` and `FindDifferentiallyExpressed(expressionData, alpha=0.05,
log2FoldChangeThreshold=1.0)` returning per-gene `DifferentialExpression` records. (A legacy
`AnalyzeDifferentialExpression` predates this Welch+BH path.)

## 1. log2 fold change (`CalculateFoldChange`)

The effect size is the log2 ratio of the two group means, **treatment over control**, with a
**pseudocount c=1** added to each mean so the ratio stays finite when a mean is 0:

```
log2FC = log2( (mean(treatment) + c) / (mean(control) + c) )      c = 1
```

- **Sign convention:** positive = **up in treatment** (condition 2 / group 2), negative =
  down (DESeq2 Love et al. 2014; Science Park RNA-seq lesson: log2FC of +1 = a 2× higher level in the
  treated state). The DOWN case is the exact negative of the UP case.
- The pseudocount only perturbs genes with a **near-zero mean**; the unregularized definition
  `log2(mean2/mean1)` is recovered as the means grow large. It never changes the sign or the relative
  ordering of non-degenerate ratios (ASSUMPTION — no source mandates a specific pseudocount value; 1
  is the standard regularization and the value the existing analyzer already uses).

## 2. Raw p-value — Welch's two-sample t-test

Per gene, an **unequal-variance (Welch) t-test** over the two groups' replicate values, using the
**unbiased (N−1) sample variance**:

```
t = (X̄₂ − X̄₁) / sqrt( s₁²/N₁ + s₂²/N₂ )          sᵢ² = Σ(x−x̄ᵢ)²/(Nᵢ−1)
ν = (s₁²/N₁ + s₂²/N₂)² / [ s₁⁴/(N₁²(N₁−1)) + s₂⁴/(N₂²(N₂−1)) ]     (Welch-Satterthwaite df)
two-sided p = I_{ν/(ν+t²)}(ν/2, 1/2)              (Student t tail via regularized incomplete beta)
```

- **No pooled variance** — the denominator handles unequal variances/sizes (Welch 1947).
- The two-sided p-value uses the exact Student-t CDF identity `F(t) = 1 − ½·I_{ν/(t²+ν)}(ν/2, ½)`
  for `t>0`, so `P(|T| ≥ t) = I_{ν/(ν+t²)}(ν/2, ½)` — computed with the **regularized incomplete beta
  function** (not a normal approximation).

## 3. FDR adjustment — Benjamini-Hochberg

Raw p-values across all `m` genes are adjusted with the **BH step-up procedure** (Benjamini &
Hochberg 1995), reproducing R `p.adjust(method="BH")` exactly:

```
sort p descending, multiply each by m/rank (rank = m..1), running minimum, clamp to 1, restore order
   pmin(1, cummin(n/i * p[o]))[ro]
```

The adjusted p-values are **monotone non-decreasing in p-order** and each is **≥ its raw p-value and
≤ 1** (BH invariant). Correcting the thousands of simultaneous per-gene tests is essential — reporting
raw p-values inflates false positives.

## 4. DE decision — the two-criterion gate

A gene is significant **only when BOTH**:

```
|log2FC| ≥ log2FoldChangeThreshold    (default 1.0)   AND   adjusted p < alpha   (default 0.05)
```

Failing **either** criterion → **not significant** (Science Park lesson; DESeq2). The alpha is the
**adjusted-p (FDR)** threshold, strict `<`.

## Method contract (algorithm spec)

Two static entry points on `TranscriptomeAnalyzer` (`Seqeron.Genomics.Annotation`, the Annotation
server), per the TRANS-DIFF-001 spec:

- `CalculateFoldChange(expression1, expression2)` → `double` log2FC (see §1). `null`/empty list ⇒
  treated as mean 0. **O(s)** time / O(1) space (two mean passes).
- `FindDifferentiallyExpressed(genes, alpha=0.05, log2FoldChangeThreshold=1.0)` →
  per-gene `DifferentialExpression` records. `genes` is
  `IEnumerable<(string GeneId, IReadOnlyList<double> Condition1, IReadOnlyList<double> Condition2)>`
  where **Condition1 = control / reference** and **Condition2 = treatment / numerator**;
  `null` enumerable ⇒ empty result. Constraints: `0 < alpha ≤ 1`, `log2FoldChangeThreshold ≥ 0`.
  **O(g·s + g·log g)** time / O(g) space — the `g·log g` term is the BH sort.

**Inputs are treated as already-normalized** expression — no internal library-size / TPM normalization
is applied (that is the upstream job of [[expression-quantification]]'s `CalculateTPM` /
`QuantileNormalize`). No exceptions are thrown for degenerate input; degenerate cases resolve to the
conventions in *Invariants and edge cases* below.

**Output record fields** (one per gene):

| Field | Type | Meaning |
|-------|------|---------|
| `GeneId` | string | gene identifier |
| `Log2FoldChange` | double | `log2((m2+c)/(m1+c))`; positive = up in condition 2 |
| `PValue` | double | raw two-sided Welch t-test p-value |
| `AdjustedPValue` | double | Benjamini-Hochberg FDR-adjusted p-value |
| `IsSignificant` | bool | `|log2FC| ≥ threshold` AND `adjusted p < alpha` |
| `Regulation` | string | `"Upregulated"` / `"Downregulated"` / `"Unchanged"` by the sign of log2FC |

This unit is **purely numeric over expression vectors** — no substring search or pattern matching, so
the repository suffix tree does not apply.

## Invariants and edge cases

- **INV:** adjusted p ∈ [raw p, 1]; DE ⇔ `|log2FC| ≥ threshold AND adjP < alpha`.
- **Fewer than 2 replicates** in a group → the (N−1) sample variance is undefined → the gene cannot be
  tested → **p = 1** (not significant). Conservative convention; never affects N≥2 inputs.
- **Zero pooled SE** (both groups constant): equal means → 0/0 → **p = 1**; unequal means → t = ±∞ →
  **p = 0** (perfect separation with zero variance).
- **Zero mean expression** → the pseudocount keeps `log2FC` finite (see §1).
- **Empty input → empty output**; identical groups → log2FC 0, p 1, not significant.

Worked oracles (from [[trans-diff-001-evidence]]):
- **Fold change** (control {10,10,10} vs treatment {40,40,40}): `log2(41/11) = +1.8981204…` (UP); the
  reversed groups give the exact negative **−1.8981204…** (DOWN); equal groups → **0** (FLAT).
- **Welch t-test** control {1,2,3} vs treatment {7,8,9}: X̄=2,8; s²=1,1; se=√(2/3)=0.8164966;
  **t = 7.348469**, **ν = 4** (equal variances/sizes), two-sided **p = 0.0018262607** — cross-checked
  against SciPy `ttest_ind([7,8,9],[1,2,3], equal_var=False)` (t=7.3484692283, p=0.0018262607), p to 1e-6.
- **BH** raw (0.001, 0.4, 0.5, 0.9) → adjusted **(0.004, 0.6667, 0.6667, 0.9)**, monotone
  non-decreasing; the degenerate raw (0.01,0.02,0.03,0.04) → all **0.04**.

## Relationship to the other two-group tests

DE is one of several **per-feature two-group significance** units in the wiki, distinguished by its
**statistical machinery**:

- **[[significant-taxa-detection]]** (metagenomics differential abundance) uses the **non-parametric
  Mann-Whitney U / Wilcoxon rank-sum** test and applies **no built-in FDR** (BH is the caller's job) —
  DE instead uses a **parametric Welch t-test on log2FC effect sizes** and **does apply BH** across
  genes. The `alternative_to` counterpart.
- **[[differentially-methylated-regions]]** (epigenetics) compares two samples with **Fisher's exact
  test** on pooled methylated/unmethylated counts, gated by a q-value + %methylation-difference cutoff
  — structurally the same "effect-size threshold AND adjusted-p threshold" two-criterion gate, but a
  different test and a different quantity (methylation, not expression).
- **[[expression-outlier-zscore-signature-score]]** (oncology) is a **single-sample** per-gene
  z-score / signature score — a within-sample outlier call, not a between-group comparison.

DE gene lists are the natural **input to over-representation / enrichment analysis**
([[pathway-enrichment-ora]]). Upstream, the two-group test operates on **normalized** expression —
the within-sample TPM/FPKM and cross-sample quantile-normalization corrections of the sibling
[[expression-quantification]] unit (TRANS-EXPR-001).

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the log2FC + Welch-t + BH
pipeline. It is a **simple two-group** estimator — **not** the full DESeq2/edgeR negative-binomial GLM
with dispersion shrinkage (DESeq2 is cited only for the log2FC definition, sign convention, and the
BH-as-standard-FDR fact, not as the fitted model). Three source-backed assumptions: pseudocount = 1;
<2 replicates → p = 1; se = 0 → p = 1 (equal means) / p = 0 (unequal means). All formulas match their
primary sources (DESeq2 log2FC/BH; Welch 1947 t-statistic; Student-t CDF regularized-incomplete-beta
identity; R `p.adjust` BH) — no source contradictions. **Not for clinical use.** The Transcriptome sibling
[[expression-quantification]] (TPM/FPKM + quantile normalization, TRANS-EXPR-001) supplies the
upstream normalized input, and [[alternative-splicing-psi]] (PSI / ΔPSI, TRANS-SPLICE-001) is the
splicing-level counterpart — a gene can be differentially *spliced* without changing its total
expression; further siblings (PCA/clustering) would enrich this family anchor.
