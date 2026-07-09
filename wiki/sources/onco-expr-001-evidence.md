---
type: source
title: "Evidence: ONCO-EXPR-001 (tumor gene-expression z-score outlier + combined-z signature score)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-EXPR-001-Evidence.md
sources:
  - docs/Evidence/ONCO-EXPR-001-Evidence.md
source_commit: 5f2e01c34bbe92d69b44839590da6713907640d3
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-EXPR-001

The validation-evidence artifact for test unit **ONCO-EXPR-001** — **Tumor Gene Expression
Outlier (z-score) and Signature Score**. The **thirteenth ingested unit of the Oncology
family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. It is the **first expression /
transcriptome method** to enter the wiki; the distinct method is synthesized in
[[expression-outlier-zscore-signature-score]]. [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the z-score formula corroborated four ways, no contradictions):**
  - **cBioPortal z-score normalization spec** (rank 5, platform spec) — per-gene z = **`(r − μ)/σ`**
    where `r` is the raw expression value and `μ`, `σ` are the mean/SD of a **reference (base)
    population**; two base-population choices: **diploid** (samples with gene copy-number = 0) or
    **all samples**. **Zero-variance → NA** (σ = 0); diploid method → NA when the gene has no diploid
    samples.
  - **cBioPortal `NormalizeExpressionLevels.java`** (rank 3, reference implementation) — resolves the
    prose spec's unstated details: **`avg` = arithmetic mean** `Σvᵢ/n`; **`std` = sample SD with
    divisor (n−1)**: `Σ(vᵢ−avg)²; /(n−1); sqrt`. **Zero SD → `fatalError` (throws RuntimeException)**
    rather than NA — this is the behaviour the unit adopts (throw, not NA).
  - **cBioPortal FAQ** (rank 5) — **default outlier threshold ±2**: "samples with expression z-scores
    **> 2** or **< −2** … are considered altered"; z > 2 ⇒ over-expressed, z < −2 ⇒ under-expressed;
    ±2 ≈ 95th percentile of a normal. The rule is **strict** (exactly ±2 is **not** an outlier).
  - **Lee et al. (2008)** *PLoS Comput Biol* 4(11):e1000217 (PMC2563693, rank 1) — the **combined
    z-score** (pathway/signature activity): per-gene standardize `zᵢⱼ = (gᵢⱼ − μᵢ)/σᵢ` (mean 0, SD 1
    over samples), then for a gene set of **k** members the activity **a = (Σᵢ zᵢ)/√k** ("√ of the
    number of member genes … to stabilize the variance of the mean").
  - **GSVA vignette** (Hänzelmann 2013, rank corroboration) — confirms the "combined z-score"
    (`zscore`) method is attributed to Lee et al. (2008): per-gene z then Σz/√k over the set.

- **Documented corner cases / failure modes:**
  - **Zero-SD reference** (constant cohort) → **throws** (per `getZ` `fatalError`); no defined z.
  - **Empty / single-sample reference** (n ≤ 1) → the (n−1) divisor is undefined; sample SD needs ≥ 2
    reference values.
  - **Threshold boundary** — strict `> 2` / `< −2`; z = ±2 exactly is **not** an outlier.
  - **Single-gene signature** (k = 1) → a = z₁/√1 = z₁ (well defined). **Empty signature** (k = 0) →
    Σz/√0 undefined → invalid input.

- **Datasets (deterministic, hand-derived):**
  - **Reference cohort (gene G) {2, 2, 4, 6, 6}:** μ = 4, Σ(rᵢ−μ)² = 16, sample variance 16/(5−1) = 4,
    **σ = 2**. Then z(x=10) = **3.0** (outlier, over), z(x=8) = **2.0** (boundary, **NOT** outlier),
    z(x=4) = **0.0**, z(x=−1) = **−2.5** (outlier, under).
  - **Signature (combined z) z = {3, 1, −1, 1}:** k = 4, Σz = 4, √k = 2, **a = 2.0**. Single-gene set
    z = {2.5} → **a = 2.5**.

- **Coverage recommendations (6 items):** MUST — per-gene z = (x−μ)/σ with **σ = sample SD (n−1)** on
  the {2,2,4,6,6} cohort; strict ±2 outlier classification (x=8 z=2.0 and x=4 NOT outliers);
  combined signature a = Σz/√k = 2.0 (and 2.5 single-gene); **zero-SD reference throws**. SHOULD —
  null/empty cohort, size-1 reference (SD undefined), empty signature (argument validation). COULD —
  symmetry/monotonicity property (z increases with x, sign flips).

## Deviations and assumptions

- **ASSUMPTION — caller supplies the reference cohort and signature gene set.** Reference distributions
  and signatures are caller inputs; the algorithm bundles no specific cohort or signature. An API/scope
  decision, not a numeric one.
- **ASSUMPTION — inputs are already on a scale where a z-score is meaningful** (e.g. log-transformed
  TPM/FPKM/microarray intensity), matching cBioPortal's "raw expression value `r`". The z formula is
  scale-agnostic; this changes no computed output for given inputs.
- **DEVIATION (behavioural) — zero-SD reference throws** (following `NormalizeExpressionLevels.java`
  `fatalError`) rather than emitting **NA** as the prose spec's "Z-Score ← NA when SD = 0" would. The
  two cBioPortal sources themselves differ on this point; the unit follows the reference
  implementation.

No further source contradictions — the z formula is corroborated four independent ways; the sole
NA-vs-throw divergence is resolved in favour of the reference code.
</content>
</invoke>
