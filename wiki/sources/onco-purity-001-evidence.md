---
type: source
title: "Evidence: ONCO-PURITY-001 (tumor purity from somatic SNV VAF — CNAqc expected-VAF inversion)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-PURITY-001-Evidence.md
sources:
  - docs/Evidence/ONCO-PURITY-001-Evidence.md
source_commit: fdf583e25989b1d2bcbc999fa056fb16119f8c31
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-PURITY-001

The validation-evidence artifact for test unit **ONCO-PURITY-001** — **Tumor Purity Estimation** from
somatic-SNV VAF / allele-specific copy number (`EstimatePurityFromVAF`, the copy-neutral diploid
heterozygous `π = 2·VAF` closed form, and `EstimatePurity`, the general inversion of the CNAqc
expected-VAF formula). The **twenty-eighth ingested unit of the Oncology family** and one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct
method is synthesized in its own concept, [[tumor-purity-from-mutation-vaf]]; [[test-unit-registry]]
tracks the unit.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **CNAqc — vignette + Genome Biology 2024 paper (Antonello et al., doi:10.1186/s13059-024-03170-5)**
    (paper rank 1, package docs rank 3) — the **expected-VAF formula (verbatim)**
    `v_m(c) = m·π·c / [2(1−π) + π(n_A + n_B)]`, with `m` = multiplicity ("mutations present in m copies
    of the tumour genome"), `π` = tumour purity, `c` = clonality / cancer-cell fraction (0<c<1), and a
    segment written `n_A:n_B` for the allele-specific copy numbers. Denominator = mean allele-copy count
    per cell: `2(1−π)` = healthy diploid normal, `π·n_tot` = tumour reads. Worked numeric examples:
    real purity 60% → expected VAF 30% for a diploid het mutation (⇒ **v = π/2, π = 2·VAF**), tolerance
    band 55–65% → VAF 27.5–32.5%; and a **2:1 segment shows two clonal peaks at 33% and 66%** (m = 1 →
    1/3, m = 2 → 2/3, n_tot = 3, π = 1). General two-state form
    `v = (m₁ρ₁ + m₂ρ₂)π / {2(1−π) + π[(n_{A,1}+n_{B,1})ρ₁ + (n_{A,2}+n_{B,2})ρ₂]}` reduces to the
    vignette form at a single clonal state (c = 1).
  - **FACETS — Shen & Seshan (NAR 2016, doi:10.1093/nar/gkw520)** (rank 1) — the mixing model
    `m* = mΦ + (1−Φ)`, `p* = pΦ + (1−Φ)` mixes a normal diploid `(1,1)` genotype with an aberrant
    `(m,p)` at cellular fraction Φ, **independently confirming the `2(1−π) + π·n_tot` denominator**
    (normal contributes 2 copies weighted 1−π).
  - **ABSOLUTE — Carter et al. (Nature Biotechnology 2012, doi:10.1038/nbt.2203)** (rank 1) — converts
    allelic fractions of point mutations into per-cancer-cell allele counts (multiplicity) "by
    correcting for sample purity and local copy-numbers" — the same purity/copy-number correction,
    inverted here to estimate purity.

- **Documented corner cases / failure modes:**
  - **Multiplicity ambiguity on amplified segments:** a 2:1 (n_tot = 3) segment yields two VAF peaks
    (1/3, 2/3) because m may be 1 or 2 — purity is not inferable from VAF alone without the CN state and
    m; copy-neutral diploid het (1:1, m = 1, n_tot = 2) avoids this, giving the robust `π = 2·VAF`.
  - **Subclonal (c < 1):** VAF depressed (v ∝ c); treating subclonal as clonal underestimates purity —
    purity must be estimated from clonal mutations.
  - **Purity < 0.1 (below detection):** VAFs near sequencing noise; high stromal contamination ≡ low
    purity; no heterozygous SNPs / informative variants → purity undefined.

- **Datasets (deterministic worked oracles):**
  - **Clonal het diploid worked example** (CNAqc): segment 1:1 (n_tot = 2), m = 1, c = 1 → expected VAF
    0.30 at purity 0.60 (= π/2); purity from VAF 0.30 = 0.60 (= 2·VAF); tolerance 0.55–0.65 ↔ VAF
    0.275–0.325.
  - **2:1 amplified segment, purity = 1** (CNAqc): n_tot = 3 → clonal VAF 1/3 (m = 1) and 2/3 (m = 2),
    both recovering π = 1 via the general inversion.

- **Coverage recommendations:** MUST test `EstimatePurityFromVAF` on VAF 0.30 → 0.60 and the boundaries
  VAF 0.50 → 1.0 / 0.0 → 0.0; MUST test `EstimatePurity` (allele-specific) recovering π = 1 on a 2:1
  segment for both m = 1 (v = 1/3) and m = 2 (v = 2/3), and agreeing with the VAF-only estimator on a
  diploid het segment (VAF 0.30 → 0.60); MUST test invalid inputs rejected (VAF outside [0,1], diploid
  VAF > 0.5 implying π > 1, empty variant list, non-positive CN); SHOULD test median aggregation across
  several mixed-VAF SNVs; COULD test purity-below-detection (VAF near 0) → small π without error.

## Deviations and assumptions

- **ASSUMPTION — VAF-only estimator uses the copy-neutral diploid heterozygous model.**
  `EstimatePurityFromVAF` assumes supplied variants are clonal (c = 1), heterozygous (m = 1) somatic
  SNVs at copy-neutral diploid (n_tot = 2) loci, giving `π = 2·VAF`. The textbook special case and the
  band CNAqc uses, stated in the API contract — a modelling-scope choice, not an invented constant; the
  formula itself is fully source-derived.
- **ASSUMPTION — cross-variant aggregation uses the median.** Multiple clonal het SNVs are combined by
  the median per-variant purity (robust to subclonal/outlier VAFs). The literature establishes the
  per-variant relation; the robust central estimator over the set is a documented,
  non-correctness-affecting aggregation policy (it does not change the single-variant formula).

No source contradictions — CNAqc (the expected-VAF formula + worked peaks), FACETS (the mixing-model
denominator), and ABSOLUTE (the inverse purity/copy-number correction) cover disjoint aspects and agree.
