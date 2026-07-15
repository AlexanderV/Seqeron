---
type: source
title: "Evidence: ONCO-IMMUNE-001 (Immune infiltration estimation — ν-SVR deconvolution + ESTIMATE ssGSEA)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-IMMUNE-001-Evidence.md
sources:
  - docs/Evidence/ONCO-IMMUNE-001-Evidence.md
source_commit: a197fb86ceeffb8de5c09005d269f020e46584f5
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-IMMUNE-001

The validation-evidence artifact for test unit **ONCO-IMMUNE-001** — **Immune Infiltration
Estimation** (tumor immune microenvironment quantification by expression deconvolution + enrichment
scoring). The **twentieth ingested unit of the Oncology family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in its own concept, [[immune-infiltration-deconvolution]]; [[test-unit-registry]] tracks
the unit.

## What this file records

- **Online sources (mutually consistent, complementary methods, no contradictions):**
  - **Newman et al. (2015)**, Nature Methods 12(5):453–457 — **CIBERSORT**: linear mixture `m = S·f`,
    **LM22** signature matrix (547 genes × 22 hematopoietic cell types), solved by **ν-SVR** (or NNLS),
    Monte-Carlo permutation p-value; validated vs flow cytometry (rank 1).
  - **Chen et al. (2018)** CIBERSORT protocol, Methods Mol Biol 1711:243–259 — the ν sweep
    **{0.25, 0.5, 0.75}**, selection by **lowest RMSE** between m and S·f, zero-clip + normalize-to-1
    output, z-score standardisation before regression, LM22 = 547 × 22 (rank 1).
  - **Yoshihara et al. (2013)**, Nature Communications 4:2612 — **ESTIMATE**: ssGSEA over 141-gene
    immune + 141-gene stromal signatures → immune / stromal / ESTIMATE scores; **tumor purity =
    cos(0.6049872018 + 0.0001467884 × ESTIMATE_score)** (rank 1). Re-retrieval (2026-06-24) confirmed
    the coefficients against three independent sources + the `tidyestimate` reference R, and the two
    domain rules: **Affymetrix-only** calibration (fit by NLS vs ABSOLUTE on TCGA Affymetrix) and
    **negative cosine → NA** (`double.NaN`).
  - **Becht et al. (2016)**, Genome Biology 17:218 — **MCP-counter**: geometric mean of cell-type
    marker genes → abundance scores (not proportions), comparable only within a cell type (rank 1).
  - **Schölkopf, Smola, Williamson & Bartlett (2000)**, Neural Computation 12(5):1207–1245 + Smola &
    Schölkopf (2004) tutorial eqs (60)–(62) — the **ν-SVR** primal/dual, primal recovery
    `f(x)=Σ(αᵢ−αᵢ*)k(xᵢ,x)+b`, and Theorem 9 (ν bounds errors / support vectors) (rank 1).
  - **CIBERSORT / LM22 licence** (Stanford) — verbatim: **no redistribution, no modification,
    non-commercial only**; LM22 gated behind registration → **not redistributable**, not embedded
    (rank 1, governs data handling).
  - **Monaco et al. (2019)**, Cell Reports 26(6):1627–1640.e7 (PMC6367568) — **ABIS-Seq** signature
    matrix, **CC BY 4.0** (permissive-with-attribution): 1296 genes × 17 immune cell types,
    mRNA-abundance-normalised (rank 1 paper + rank 5 dataset).
  - Scikit-learn 1.6.1 `NuSVR` (libsvm) used as the cross-implementation ν-SVR reference (rank —
    tooling).

- **Implemented surface:** `DeconvoluteImmuneCellsNuSvr` (opt-in CIBERSORT ν-SVR: sweep ν, lowest-RMSE,
  zero-clip + normalize), `DeconvoluteImmuneCells` (NNLS/LLSR baseline, Abbas 2009, unchanged),
  `EstimateInfiltration` (**simplified ssGSEA** = rank-weighted mean of signature-gene expression),
  `EstimateTumorPurity` (opt-in Affymetrix-only cosine transform, negative → NaN), `LoadSignatureMatrix`
  (LM22-format TSV loader, caller-supplied), `LoadBundledAbisSignatureMatrix` (bundled ABIS-Seq).

- **Documented corner cases / failure modes:** low signal-to-noise degrades accuracy; missing
  signature genes → unreliable; **collinear cell types** (resting vs activated NK) recover support but
  looser proportions; non-hematopoietic content lands in the residual; empty expression profile → no
  scores; no overlapping genes → zero/undefined; negative (log) expression is valid; ESTIMATE score may
  exceed the valid purity range at extreme low purity; ESTIMATE inappropriate for hematological
  (non-solid) tumors.

- **Datasets (deterministic worked oracles):**
  - **Dataset 4 (ν-SVR planted truth, `DefaultSignatureMatrix`):** planted f = {CD8 0.60, B_naive 0.30,
    Monocytes 0.10} → recovered {0.5971, 0.2989, 0.1040}, others ≈ 0, Σ = 1, corr 0.99997.
  - **Dataset 5 (scikit-learn cross-check):** 3-type standardized problem → selected **ν = 0.75**,
    fractions {A 0.5085, B 0.1795, C 0.3120} vs sklearn {0.5085, 0.1795, 0.3120} (< 2×10⁻³).
  - **Dataset 6 (bundled ABIS-Seq, 1296 × 17):** planted {NK 0.60, Monocytes C 0.40} → {≈0.650,
    ≈0.350}, 15 absent types exactly 0, corr ≈ 0.996 (tol 0.06); pure {Monocytes C 1.0} → 1.0 exactly.
    Reference values: `S1PR3`/Monocytes C = 45.720735005602499; `CD8A`/T CD8 Memory = 1060.1507652944399;
    `MS4A1`/B Naive = 3220.5650656491198.
  - Datasets 1–3: pure cell type → 100% that type; equal 50:50 mixture → 0.5 each; ESTIMATE ordering
    (high-immune sample score > low-immune sample score).

- **Coverage recommendations:** MUST — empty profile → empty; single type → 100%; equal mixture →
  ~50:50; Σf = 1 and f ≥ 0 invariants; immune-score ordering; no-overlap → zero; ν-SVR planted-truth
  recovery; scikit-learn `NuSVR` cross-match; selected ν ∈ {0.25,0.5,0.75} by lowest RMSE; LM22-format
  loader parse + reject malformed; ABIS bundled dimensions + exact values + planted-truth recovery.

## Deviations and assumptions

- **LM22 caller-supplied (licence):** Stanford forbids LM22 redistribution → not embedded; only the
  ν-SVR algorithm + LM22-format loader ship. Exact CIBERSORT-LM22 parity is **not claimed** (also needs
  the tool's quantile-normalisation pipeline). *Partially resolved* by bundling the permissive ABIS-Seq
  matrix (Monaco 2019, CC BY 4.0) so deconvolution works out-of-the-box — but ABIS (17 types) is not an
  LM22 substitute.
- **ν-SVR resolved:** `DeconvoluteImmuneCellsNuSvr` implements the CIBERSORT ν-SVR (verified by
  planted-truth + scikit-learn `NuSVR` cross-check); the NNLS/LLSR baseline is retained unchanged.
- **Simplified ssGSEA (assumption):** `EstimateInfiltration` uses a rank-weighted mean of signature
  genes rather than the full GSEA enrichment-score computation — retains the essential scoring concept.
- **ESTIMATE purity domain:** the cosine transform is opt-in, **Affymetrix-only** (per the reference
  implementation's `is_affymetrix` gate), and negative cosine values → **NaN** (out of calibrated
  domain). cos is monotone-decreasing on [0, π] → purity strictly decreases with ESTIMATE score.

No source contradictions — CIBERSORT (ν-SVR deconvolution), ESTIMATE (ssGSEA scoring), MCP-counter
(marker geometric mean), and the Schölkopf ν-SVR formulation cover complementary methods and agree on
the shared linear-mixture / signature framing.
