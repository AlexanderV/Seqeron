---
type: concept
title: "Differentially methylated region (DMR) detection"
tags: [epigenetics, algorithm]
sources:
  - docs/Evidence/EPIGEN-DMR-001-Evidence.md
  - docs/algorithms/Epigenetics/Differentially_Methylated_Regions.md
source_commit: 547b4cdb6c5e5e4454d7cbdd9005a4741adc00e0
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: epigen-dmr-001-evidence
      evidence: "Test Unit ID: EPIGEN-DMR-001 ... Algorithm: Differentially Methylated Region (DMR) detection (tiling-window + Fisher's exact test, methylKit model)"
      confidence: high
      status: current
---

# Differentially methylated region (DMR) detection

Finding **genomic regions whose DNA methylation differs between two samples** — the methylKit
(Akalin et al. 2012) **tiling-window + Fisher's-exact-test** model. This is the **fifth ingested
unit of the Epigenetics family** and a **distinct algorithm** from its siblings: it *compares*
methylation between two groups rather than measuring a single sample. It **consumes the per-CpG
β-values** that [[bisulfite-methylation-calling]] produces (methylation = C/(C+T) per base), over
the CpG sites that [[cpg-island-detection]] locates; contrast the single-sample
[[epigenetic-age-horvath-clock]] (age from β-values) and the histone-mark
[[chromatin-state-prediction]]. Validated under test unit **EPIGEN-DMR-001**; the record is
[[epigen-dmr-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

A **DMR** is "a base/region ... classified as a differentially methylated cytosine (DMC) or region
(DMR)" once the null hypothesis of equal methylation is rejected — i.e. a region of **adjacent CpG
sites** that are differentially methylated between the two samples (Akalin 2012).

## 1. Tiling windows (`tileMethylCounts`)

The genome is partitioned into fixed **tiling windows** and per-base methylated/unmethylated counts
are summed within each window. methylKit defaults: **`win.size=1000`, `step.size=1000`,
`cov.bases=0`** (non-overlapping 1 kb tiles; `cov.bases` filters tiles with too few covered bases,
default 0 keeps all). Sites farther apart than `win.size` fall into **separate** windows; a tile must
contain **≥1 informative cytosine** to be tested.

## 2. Per-window methylation difference

Per-site methylation is the Bismark-style ratio **C/(C+T)** (methylated reads over total covered
reads). For each window the signed difference is

```
meth.diff = group2 %methylation − group1 %methylation      (percentage points)
```

## 3. Significance — Fisher's exact test (one sample per group)

With **one sample per group** (the implemented regime — e.g. after pooling), methylKit applies
**Fisher's exact test**; with replicates it uses logistic regression instead (out of scope here). The
covered cytosines inside a tile are **pooled into one 2×2 table** — summed methylated vs unmethylated
reads per group — and Fisher's exact test is applied to the window. The single-table hypergeometric
probability for cells a,b,c,d (row sums a+b, c+d; column sums a+c, b+d; total n) is

```
p = (a+b)! (c+d)! (a+c)! (b+d)! / (a! b! c! d! n!)
```

The **two-sided p-value** sums the probabilities of all tables (same fixed margins) whose
hypergeometric probability ≤ that of the observed table. A q-value adjusts for multiple testing.

## 4. Selection and hyper/hypo classification (`getMethylDiff`)

Default reporting extracts regions with **q-value < 0.01 AND |meth.diff| > 25** (`get.methylDiff`
signature `difference=25, qvalue=0.01, type="all"`; `difference` is in **percentage points**). The
threshold is **strict greater-than**:

- **Hyper-methylated** ← `qvalue < 0.01 & meth.diff > +difference` (higher methylation than control).
- **Hypo-methylated** ← `qvalue < 0.01 & meth.diff < −difference` (lower methylation than control).

A window whose **|meth.diff| is at or under the cutoff is NOT reported** (25 exactly → excluded).

## Invariants and edge cases

- **INV:** `meth.diff = group2% − group1%` ∈ [−100, +100]; hyper ⇔ meth.diff > 25, hypo ⇔ < −25,
  both gated by q < 0.01 (strict `>`/`<`).
- **Empty input** (no covered positions) → no tiles, **no DMRs** (empty result set).
- **Zero-coverage group** in a window (a row total of 0) → the 2×2 table is degenerate, the only
  feasible table given the margins is the observed one, so the exact two-sided **p = 1.0** (not
  reported).
- **Identical proportions / zero marginal** (a row or column total 0) → only one feasible table →
  **p = 1.0**.
- **Adjacency:** a DMR is a region of *adjacent* CpG sites — a window with too few covered sites
  (`cov.bases` / `minCpGCount`) is not reported.
- **Determinism:** same inputs → identical DMR list ordered by start position.

Worked oracles: the Fisher single-table probability for a=1, b=9, c=11, d=3 (n=24) is
**≈0.001346076** (Wikipedia "studying by gender" example, validating the hypergeometric term);
a hyper window with group1 fully unmethylated (level 0.0) and group2 fully methylated (level 1.0),
coverage 20 across 3 sites, pools to methylated {g1=0, g2=60} / unmethylated {g1=60, g2=0} →
**meth.diff = +100**, two-sided Fisher **p ≈ 0** (complete separation), classified
**Hypermethylated**.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the tile→pool→Fisher
pipeline, **not** a full methylKit port. **Not implemented:** the logistic-regression (replicate)
path; SLIM/other q-value models beyond BH-style adjustment; overlapping-window (`step<win`) tiling.
Two evidence-backed **assumptions**: (1) per-window pooling of single-sample sites into one 2×2 table
mirrors `tileMethylCounts` (sum counts over the window) followed by Fisher's exact on the tiled
counts — the documented methylKit tile→test pipeline; (2) integer numC/numT are reconstructed as
`round(level × coverage)` and `round((1−level) × coverage)` because the repository `MethylationSite`
stores a fractional level + integer coverage rather than raw C/T counts (a representation detail, not
a change to the C/(C+T) definition). No source contradictions — Akalin 2012, the methylKit
`tileMethylCounts`/`calculateDiffMeth`/`getMethylDiff` reference, and the Fisher's-exact-test
hypergeometric definition are mutually consistent.
