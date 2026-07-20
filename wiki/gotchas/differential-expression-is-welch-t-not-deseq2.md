---
type: gotcha
title: "differential_expression is a two-group Welch t-test, not DESeq2/edgeR"
tags: [transcriptome, gotcha]
mcp_tools:
  - differential_expression
sources:
  - docs/algorithms/Transcriptome/Differential_Expression.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# differential_expression is a two-group Welch t-test, not DESeq2/edgeR

**The trap.** `differential_expression` looks like an RNA-seq DE caller, but it is a **simple
two-group estimator**: per-gene log2 fold change (pseudocount 1) + **Welch's t-test** + a
Benjamini–Hochberg FDR pass. It is **not** the DESeq2/edgeR negative-binomial GLM with
dispersion shrinkage — DESeq2 is cited only for the log2FC definition, the sign convention, and
the BH-as-standard-FDR fact, **not** as the fitted model.

**Why it bites.** With the few replicates typical of real RNA-seq (n = 2–3), a Welch t-test has
**no dispersion shrinkage across genes**, so small-N significance is unreliable and will not match
DESeq2/edgeR on the same counts. Three conservative degeneracies also surprise a caller:

- **< 2 replicates in a group** → the (N−1) sample variance is undefined → the gene is
  **untestable → p = 1** (not "not expressed").
- **Zero pooled SE** (both groups constant): equal means → 0/0 → **p = 1**; unequal means →
  t = ±∞ → **p = 0** (perfect separation reported as maximally significant).
- The input must already be **normalized** — there is no library-size / size-factor model here.
  [[expression-quantification]] (TPM/FPKM + quantile normalization) supplies that upstream input.

**What to rely on instead.** For count-based small-N RNA-seq that needs DESeq2-grade calls, fit a
real negative-binomial GLM externally. Use this tool for a fast, validated two-group screen on
**already-normalized** values, and read the BH-adjusted p as a screen, not a clinical verdict
(research-grade, [[research-grade-limitations]]). Full model: [[differential-expression]]; the
splicing-level counterpart is [[alternative-splicing-psi]].
