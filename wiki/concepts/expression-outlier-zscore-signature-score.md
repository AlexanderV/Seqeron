---
type: concept
title: "Gene-expression z-score outlier + combined-z signature score"
tags: [oncology, transcriptome, algorithm]
sources:
  - docs/Evidence/ONCO-EXPR-001-Evidence.md
  - docs/Evidence/ONCO-IMMUNE-001-Evidence.md
source_commit: a197fb86ceeffb8de5c09005d269f020e46584f5
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-expr-001-evidence
      evidence: "Test Unit ID: ONCO-EXPR-001, Algorithm: Tumor Gene Expression Outlier (z-score) and Signature Score"
      confidence: high
      status: current
---

# Gene-expression z-score outlier + combined-z signature score

The **thirteenth ingested Oncology unit** (ONCO-EXPR-001) and the wiki's **first
expression / transcriptome method**. It turns a tumor's per-gene expression into two
interpretable quantities: a **per-gene outlier z-score** (is this gene over- or
under-expressed relative to a reference cohort?) and a **combined-z signature/pathway
activity score** (how active is a gene set as a whole?). Validated under test unit
**ONCO-EXPR-001** ([[onco-expr-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern. The numerical core is
**exact** with respect to the cited formulas.

Both quantities are the standard **cBioPortal** (outlier) and **Lee et al. 2008 / GSVA**
(signature) statistics — see [[onco-expr-001-evidence]] for the source-by-source trace.

## Per-gene outlier z-score (cBioPortal)

For a query expression value `r` against a **reference (base) population**:

```
z = (r − μ) / σ
```

- **μ** = arithmetic mean of the reference values; **σ** = **sample standard deviation with
  divisor (n − 1)** (`Σ(vᵢ−μ)²/(n−1)`, then √). The reference-implementation `std()` settles
  the n-vs-(n−1) question the prose spec leaves open — it is **(n−1)**.
- **Reference population** is caller-supplied and is one of cBioPortal's two conventions:
  the **diploid** samples (gene copy-number = 0) or **all** samples with expression values.
- **Outlier classification** uses the default **±2** threshold, **strict**: `z > 2` →
  over-expressed, `z < −2` → under-expressed, and **exactly ±2 is NOT an outlier**. ±2 ≈ the
  95th percentile of a normal (two SDs from the mean).

## Combined-z signature / pathway activity (Lee et al. 2008)

Standardize each member gene to a z-score, then average over the **k** genes of the set with
a **√k** (not `k`) denominator to stabilize the variance of the mean:

```
a = (Σᵢ zᵢ) / √k
```

This is GSVA's "combined z-score" (`zscore`) method, attributed to Lee et al. (2008). A
single-gene set (k = 1) reduces to `a = z₁`.

## Invariants and edge cases

- **Zero-SD reference throws.** A constant (degenerate) reference cohort has σ = 0 and no
  defined z — the reference implementation `fatalError`s, so the unit **throws** (a
  behavioural deviation from cBioPortal's prose "z ← NA when SD = 0"; the two cBioPortal
  sources disagree, and the code wins — see [[onco-expr-001-evidence]]).
- **Reference size ≥ 2.** The (n−1) divisor is undefined for n ≤ 1; a single reference value
  has no spread.
- **Empty signature (k = 0)** → Σz/√0 undefined → invalid input; **k = 1** is well defined.
- **Scale assumption.** Inputs are assumed already on a scale where a z is meaningful
  (log-transformed TPM/FPKM/microarray intensity); the formula itself is scale-agnostic.

Worked oracles (from [[onco-expr-001-evidence]]): reference `{2,2,4,6,6}` → μ = 4, σ = 2, so
`x=10 → 3.0` (over), `x=8 → 2.0` (boundary, not outlier), `x=4 → 0.0`, `x=−1 → −2.5` (under);
signature `z = {3,1,−1,1}` → `a = 4/2 = 2.0`, single-gene `{2.5}` → `a = 2.5`.

## Scope and limitations

A [[research-grade-limitations|research-grade]] statistic, **not for clinical use**. It
computes z and activity from **caller-supplied** normalized expression, references, and gene
sets — it bundles no cohort, signature, or gene-set database, and applies no
multiple-testing correction across genes. It is the **expression-layer** Oncology unit,
orthogonal to the copy-number / clonal-structure ONCO units (e.g.
[[copy-number-alteration-classification]], [[cancer-cell-fraction-clonal-clustering]]); the
diploid reference-population option intersects the copy-number world only as an input filter.
Unlike the ranking-based GSEA, the combined-z score treats all set members equally and the
outlier test is threshold-based. It shares the **single-sample signature/pathway-activity scoring**
layer with the ssGSEA immune/stromal scores of [[immune-infiltration-deconvolution]] (ONCO-IMMUNE-001,
the other expression/transcriptome ONCO unit). No source contradictions beyond the single NA-vs-throw
divergence, resolved to throw.
