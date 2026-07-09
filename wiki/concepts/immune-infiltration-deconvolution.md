---
type: concept
title: "Immune infiltration estimation — expression deconvolution (ν-SVR / ESTIMATE)"
tags: [oncology, transcriptome, algorithm]
sources:
  - docs/Evidence/ONCO-IMMUNE-001-Evidence.md
source_commit: a197fb86ceeffb8de5c09005d269f020e46584f5
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-immune-001-evidence
      evidence: "Test Unit ID: ONCO-IMMUNE-001, Algorithm: Immune Infiltration Estimation"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:expression-outlier-zscore-signature-score
      source: onco-immune-001-evidence
      evidence: "ESTIMATE uses single-sample GSEA (ssGSEA) to compute enrichment scores for immune and stromal gene signatures — the same signature/pathway-activity scoring layer over normalized bulk expression as the combined-z signature score."
      confidence: medium
      status: current
---

# Immune infiltration estimation — expression deconvolution (ν-SVR / ESTIMATE)

The **twentieth ingested Oncology unit** (**ONCO-IMMUNE-001**) and the second wiki
**expression / transcriptome** method (after [[expression-outlier-zscore-signature-score]]).
It quantifies the **tumor immune microenvironment** from bulk RNA expression along two
complementary axes:

1. **Cell-type deconvolution** — estimate the *fractions* of immune cell types that make up a
   bulk mixture (CIBERSORT-style **ν-SVR**, plus an NNLS/LLSR baseline).
2. **Immune / stromal enrichment scoring** — a single *score* for how immune- or stromal-rich a
   sample is (**ESTIMATE** ssGSEA), and a derived **tumor purity** estimate.

The literature-traced record is [[onco-immune-001-evidence]]; [[test-unit-registry]] tracks the
unit and [[algorithm-validation-evidence]] describes the evidence-artifact pattern. This is the
antigen-agnostic quantitative sibling of the HLA / neoantigen-presentation unit
[[hla-nomenclature-and-allele-specific-loh]] in the immuno-oncology area.

## The linear mixture model (deconvolution)

Bulk expression is modelled as a non-negative mixture of reference cell-type signatures:

```
m = S · f
```

where **m** is the mixture expression vector (one value per gene), **S** the **signature matrix**
(genes × cell types), and **f** the cell-type **fraction** vector. Fractions are constrained
**f ≥ 0** and normalized to **Σf = 1**.

### ν-SVR deconvolution — `DeconvoluteImmuneCellsNuSvr` (CIBERSORT)

The CIBERSORT method (Newman et al. 2015) solves the mixture with **ν-support-vector regression**
(Schölkopf et al. 2000; linear kernel), which is robust to noise and unknown/non-immune content:

- **z-score standardize** both the mixture **m** and each signature column of **S** before regression
  (Newman 2015 Online Methods).
- **Sweep ν ∈ {0.25, 0.5, 0.75}** and select the ν giving the **lowest RMSE** between the observed
  mixture **m** and the reconstruction **S·f** (Chen et al. 2018 protocol).
- **Zero-clip** the regression coefficients (negatives → 0) and **normalize** the survivors to
  sum 1 → the reported fractions.

The engine was cross-validated two ways: **planted-truth recovery** (`m = S·f` for known f) and a
**cross-implementation match against scikit-learn 1.6.1 `NuSVR`** (libsvm backend). ν's semantics
(Schölkopf Theorem 9): ν upper-bounds the fraction of errors and lower-bounds the fraction of
support vectors.

### NNLS / LLSR baseline — `DeconvoluteImmuneCells`

The default (non-ν-SVR) path retains the older non-negative-least-squares / linear-least-squares
regression baseline (Abbas et al. 2009): `min ‖m − S·f‖²` s.t. `f ≥ 0`, then normalize to Σf = 1.
CIBERSORT's ν-SVR outperformed NNLS/QP/LLSR in the original benchmark against flow-cytometry ground
truth. A **Monte-Carlo permutation p-value** accompanies each CIBERSORT deconvolution.

## Signature matrices (the S in m = S·f)

| Matrix | Dimensions | Provenance | Bundled? |
|--------|-----------|------------|----------|
| **LM22** (CIBERSORT) | 547 genes × **22** hematopoietic cell types | Newman 2015 | **No** — caller-supplied |
| **ABIS-Seq** | 1296 genes × **17** immune cell types | Monaco 2019, **CC BY 4.0** | **Yes** — embedded resource |
| `DefaultSignatureMatrix` | 5 markers × 22 cell types | synthetic | Yes (tests / non-LM22 default) |

**LM22 is NOT redistributable** — the Stanford CIBERSORT licence forbids redistribution,
modification, and commercial/for-profit use, and gates `LM22.txt` behind registration. So the
library implements the **ν-SVR algorithm** and an **LM22-format loader** (`LoadSignatureMatrix`)
but does **not embed** LM22; the caller supplies `LM22.txt` under their own licence. Exact
CIBERSORT-LM22 parity is **not claimed** (it also needs the tool's full quantile-normalisation
pipeline). The permissively-licensed **ABIS-Seq** matrix (Monaco et al. 2019) *is* bundled as
`Resources/ABIS_sigmatrixRNAseq.tsv` and exposed via `LoadBundledAbisSignatureMatrix()`, so
`DeconvoluteImmuneCellsNuSvr` works out-of-the-box — but ABIS (17 types, mRNA-abundance scale) is
**not** an LM22 substitute.

## ESTIMATE enrichment scoring — `EstimateInfiltration` / `EstimateTumorPurity`

ESTIMATE (Yoshihara et al. 2013) scores stromal / immune content with **single-sample GSEA
(ssGSEA)** over two 141-gene signatures, producing an **immune score**, **stromal score**, and
**ESTIMATE score = immune + stromal**. The library uses a **simplified ssGSEA** (rank-weighted mean
of signature-gene expression) — an intentional deviation that keeps the ordering semantics but not
the full Kolmogorov–Smirnov running-sum enrichment statistic. This is the same signature/pathway
scoring layer as the combined-z score of [[expression-outlier-zscore-signature-score]].

**Tumor purity** is derived from the ESTIMATE score by the verbatim cosine transform:

```
Tumor_purity = cos(0.6049872018 + 0.0001467884 × ESTIMATE_score)
```

The coefficients were confirmed against **three independent sources** plus the `tidyestimate`
reference R. Two source-decisive domain rules: (a) the model was fit by nonlinear least squares
against **ABSOLUTE** purity on **TCGA Affymetrix** data, so it is valid **only for Affymetrix-derived
scores** (not RNA-seq) — this is an **opt-in** absolute-purity method; (b) a **negative cosine value
is out of the calibrated domain → NA** (mirrored as `double.NaN`). Because cos is monotone-decreasing
on [0, π], purity **strictly decreases** as the ESTIMATE score increases over the calibrated range.

### MCP-counter (abundance scoring)

Becht et al. (2016) MCP-counter computes a **geometric mean** of cell-type-specific marker genes,
producing **abundance scores** (not proportions) that are only comparable *within* a cell type
across samples — a third scoring flavor noted in the evidence but distinct from both deconvolution
fractions and ESTIMATE scores.

## Worked oracles (from [[onco-immune-001-evidence]])

- **Planted-truth ν-SVR** (`DefaultSignatureMatrix`): planted f = {CD8 0.60, B_naive 0.30, Monocytes
  0.10} → recovered {0.5971, 0.2989, 0.1040}, all others ≈ 0, Σ = 1, correlation 0.99997.
- **scikit-learn cross-check**: 3-type standardized problem → selected **ν = 0.75** (lowest RMSE),
  fractions {A 0.5085, B 0.1795, C 0.3120} vs sklearn {0.5085, 0.1795, 0.3120} (agreement < 2×10⁻³).
- **Bundled ABIS-Seq** (1296 × 17): planted {NK 0.60, Monocytes C 0.40} → {≈0.650, ≈0.350}, 15 absent
  types exactly 0, correlation ≈ 0.996 (tolerance 0.06); pure {Monocytes C 1.0} → 1.0 exactly.
- **Reference values** (ABIS Table S5, exact): `S1PR3`/Monocytes C = 45.720735005602499;
  `CD8A`/T CD8 Memory = 1060.1507652944399; `MS4A1`/B Naive = 3220.5650656491198.

## Corner cases and limitations

- **Collinear cell types** (e.g. resting vs activated NK, T-cell-heavy mixtures) recover the correct
  *support* but looser *proportions* — a documented property of the unmodified ν-SVR engine (Newman
  2015 corner case), so planted-truth tests use well-separated lineages + a one-hot case.
- **Missing / non-overlapping signature genes** → unreliable or zero/undefined fractions; **empty
  expression profile** → no scores. **Low signal-to-noise** degrades accuracy.
- **Non-hematopoietic content** is *not* estimated by deconvolution — it falls into the residual.
- **ESTIMATE domain:** designed for **solid tumors**; hematological malignancies are inappropriate
  inputs, and the ESTIMATE score may exceed the valid purity range at extreme low-purity values.
- **Negative (log-transformed) expression values** are valid inputs.

## Scope and assumptions

A [[scientific-rigor|research-grade]] statistic, **not for clinical or diagnostic use**. Three
source-consistent scope facts: LM22 is caller-supplied (Stanford licence — no exact CIBERSORT
parity claimed); ABIS-Seq is bundled under CC BY 4.0 as a working default; and ESTIMATE uses a
**simplified ssGSEA** (rank-weighted mean, not the full enrichment-score computation). No source
contradictions — the CIBERSORT, ESTIMATE, MCP-counter, and ν-SVR sources cover complementary
methods and agree on the shared linear-mixture / signature-scoring framing.
</content>
</invoke>
