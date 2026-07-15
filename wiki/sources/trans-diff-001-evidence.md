---
type: source
title: "Evidence: TRANS-DIFF-001 (Differential expression — log2FC + Welch t-test + BH FDR)"
tags: [validation, transcriptome]
doc_path: docs/Evidence/TRANS-DIFF-001-Evidence.md
sources:
  - docs/Evidence/TRANS-DIFF-001-Evidence.md
source_commit: e00919fd4a7a3e5c624a134af281d66dfe6d4831
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: TRANS-DIFF-001

The validation-evidence artifact for test unit **TRANS-DIFF-001** — **two-group RNA-seq differential
expression**: log2 fold change + Welch's (unequal-variance) two-sample t-test + Benjamini-Hochberg
FDR. The **first ingested unit of the Transcriptome / RNA-seq family** and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The algorithm is
synthesized in its own concept, [[differential-expression]]; [[test-unit-registry]] tracks the unit.
Impl `TranscriptomeAnalyzer.CalculateFoldChange` / `FindDifferentiallyExpressed`
(`Seqeron.Genomics.Annotation`).

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Love, Huber & Anders (2014) — DESeq2**, *Genome Biology* 15:550 (PMC4302049, rank 1) — the fold
    change is the **log2 ratio between treatment and control** (coefficient sign = treatment vs
    control), reported on the **log2 scale**; Wald-test p-values are **adjusted for multiple testing by
    the Benjamini-Hochberg procedure** (establishes BH as the standard RNA-seq FDR).
  - **Science Park Study Group RNA-seq lesson** (rank 3) — `log2FC = log2(condition A / condition B)`
    with the treatment as numerator; positive = upregulated in treatment; a gene is DE only under
    **two simultaneous criteria** — `|log2FC| ≥ threshold` (commonly |1| or |2|) **AND** adjusted
    p-value `< alpha` (typically 0.01/0.001).
  - **Welch's t-test** — Wikipedia citing Welch 1947 (rank 4, formulas) —
    `t = (X̄₁−X̄₂)/√(s₁²/N₁ + s₂²/N₂)` with **corrected (N−1)** sample SDs and **no pooled variance**;
    Welch-Satterthwaite df `ν ≈ (s₁²/N₁+s₂²/N₂)² / [s₁⁴/(N₁²(N₁−1)) + s₂⁴/(N₂²(N₂−1))]`.
  - **Student's t-distribution CDF** — Wikipedia (rank 4) — two-sided tail
    `P(|T| ≥ t) = I_{ν/(ν+t²)}(ν/2, ½)` via the **regularized incomplete beta** function.
  - **Benjamini & Hochberg (1995)** — Wikipedia "False discovery rate" + R `p.adjust` manual (rank
    2–4) — BH step-up: with ascending `P(1)≤…≤P(m)`, find the largest k with `P(k) ≤ (k/m)α`; the R
    reference-implementation adjusted p-value is `pmin(1, cummin(n/i * p[o]))[ro]` (sort p descending,
    ×`n/rank`, running minimum, clamp to 1, restore order); `p.adjust(method="BH")` cites BH 1995.

- **Documented corner cases / failure modes:** two-criterion gate (fail either ⇒ not significant);
  many simultaneous tests must be FDR-adjusted; **<2 replicates per group** → (N−1) variance undefined
  → gene not testable → **p = 1**; **se = 0** → equal means ⇒ p = 1, unequal means ⇒ t = ∞ ⇒ p = 0;
  **zero mean expression** → pseudocount c=1 added to both means before the ratio (log2FC finite).

- **Datasets (documented oracles):**
  - *Fold change* — control {10,10,10} vs treatment {40,40,40} → `log2(41/11) = +1.8981204…` (UP);
    reversed → **−1.8981204…** (DOWN, exact negative); equal → **0** (FLAT). log2FC =
    `log2((mean2+1)/(mean1+1))`, pseudocount 1.
  - *Welch t-test* — control {1,2,3} vs treatment {7,8,9} → X̄ 2,8; s² 1,1; se √(2/3)=0.8164966;
    **t = 7.348469**, **ν = 4**, two-sided **p = 0.0018262607** (cross-checked vs SciPy
    `ttest_ind(equal_var=False)` t=7.3484692283 / p=0.0018262607; MUST asserts exact t + df, p to 1e-6,
    p < alpha).
  - *Benjamini-Hochberg* — raw (0.001, 0.4, 0.5, 0.9) → adjusted **(0.004, 0.6667, 0.6667, 0.9)**,
    monotone non-decreasing; degenerate raw (0.01,0.02,0.03,0.04) → all **0.04**.

- **Test-coverage recommendations:** MUST — `CalculateFoldChange` UP/DOWN/FLAT values + sign
  convention (positive ⇒ higher in treatment); Welch t/df/p for {1,2,3} vs {7,8,9}; BH reproduces R
  `p.adjust`; DE gate requires BOTH `|log2FC| ≥ threshold` AND adjusted p < alpha. SHOULD — empty →
  empty; identical groups → log2FC 0/p 1/not significant; <2 replicates → p 1. COULD — adjusted p ≥
  raw p and ≤ 1 (BH invariant).

## Deviations and assumptions

- **ASSUMPTION (pseudocount = 1):** `log2(mean2/mean1)` is undefined at a zero mean; adding 1 to both
  means is the standard regularization and the value already used by the analyzer — affects only
  near-zero-mean genes, never the sign/ordering; the unregularized DESeq2/Science Park definition is
  recovered as means grow.
- **ASSUMPTION (<2 replicates → p = 1):** Welch's statistic needs the unbiased (N−1) variance,
  undefined for N<2; emitting p=1 (not testable ⇒ not significant) is the conservative convention and
  matches the analyzer; never affects N≥2 inputs.
- **ASSUMPTION (se = 0):** equal means ⇒ p = 1 (0/0, no evidence); unequal means ⇒ p = 0 (±∞, perfect
  separation) — follows directly from the t-statistic limit.

No source contradictions — DESeq2 (log2FC/BH), Welch 1947, the Student-t CDF identity, and R
`p.adjust` BH are mutually consistent. This is a **simple two-group** estimator, not the full DESeq2
negative-binomial GLM. Research-grade, not for clinical use.
