---
type: source
title: "Evidence: ONCO-VAF-001 (variant allele frequency — empirical VAF, Wilson binomial CI, purity/ploidy correction)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-VAF-001-Evidence.md
sources:
  - docs/Evidence/ONCO-VAF-001-Evidence.md
source_commit: 68661290101fe6f70f2c89ccf5f5076fff5940ce
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-VAF-001

The validation-evidence artifact for test unit **ONCO-VAF-001** — **Variant Allele Frequency Analysis**:
the empirical `VAF = altReads/totalReads`, its **Wilson score binomial confidence interval**, and the
CNAqc **`AdjustVAFForPurity`** purity/ploidy correction. The **thirty-sixth ingested unit of the Oncology
family** and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern. The distinct method is synthesized in its own concept,
[[variant-allele-frequency-and-binomial-ci]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **GATK Mutect2 FAQ (rank 3)** — the `AF` FORMAT field is documented as "allele fractions of alternate
    alleles in the tumor" but is a **model estimate** (posterior marginalised over allele fractions), **not**
    the raw ratio. The **empirical** VAF is the model-free `altReads/totalReads` = alt AD / sum AD (GATK `AD`
    field / samtools mpileup) — the textbook definition this unit implements.
  - **Wikipedia — Binomial proportion confidence interval / Wilson 1927 (rank 4, citing the primary Wilson
    1927 JASA 22(158):209–212)** — the **Wilson score interval** verbatim: `center = (p̂ + z²/2n)/(1 + z²/n)`,
    `margin = (z/(1+z²/n))·√(p̂(1−p̂)/n + z²/4n²)`; **z = 1.96** for 95%; contrasted with the **Wald**
    interval `p̂ ± (z/√n)·√(p̂(1−p̂))` which is criticised for overshoot (bounds outside [0,1]) and zero-width
    at extreme proportions — both avoided by Wilson.
  - **Tarabichi et al. 2017/2021, PMC5538405 (rank 1)** — VAF/CCF/purity relationship; the diploid
    heterozygous clonal case gives **expected VAF = purity × 0.5** (80% purity → VAF 0.4).
  - **CNAqc — Househam/Caravagna, *Genome Biology* 2024 25:38, doi:10.1186/s13059-024-03170-5 (rank 1)** —
    the expected clonal VAF `v = m·π / (2(1−π) + π·n_tot)` (normal ploidy 2), with peaks at 33%/66% for a
    2:1 (n_tot=3) segment; the **inversion** `m·CCF = VAF·(2(1−π)+π·n_tot)/π` is the purity/ploidy
    adjustment computed by `AdjustVAFForPurity`.

- **Documented corner cases / failure modes:**
  - **VAF > 1 (alignment artifacts):** misaligned/overlapping reads can push alt above total; empirical VAF
    requires `altReads ≤ totalReads` — a violation is invalid input.
  - **Zero coverage (`totalReads = 0`):** VAF undefined (0/0); division-by-zero must be guarded.
  - **Extreme proportions (`p̂ = 0` or `1`):** Wilson stays in [0,1] with non-zero width (Wald overshoots
    / collapses).
  - **Purity `π = 0`:** the correction divides by π; pure-normal makes the adjusted fraction undefined.

- **Datasets (deterministic worked oracles):**
  - **Empirical VAF:** 25/100→0.25, 50/100→0.50, 5/20→0.25, 0/10→0.00, 10/10→1.00.
  - **Wilson 95% CI (z=1.96):** 25/100→center 0.2592, [0.1755, 0.3430]; 50/100→0.5, [0.4038, 0.5962];
    5/20→0.2903, [0.1119, 0.4687]; **0/10→0.1388, [0.0000, 0.2775]** (lower=0 no-overshoot);
    **10/10→0.8612, [0.7225, 1.0000]** (upper=1 no-overshoot).
  - **Purity/ploidy-corrected VAF (`m·CCF = VAF·(2(1−π)+π·n_tot)/π`):** (0.40, π 0.80, n_tot 2)→1.00,
    (0.35, 0.70, 2)→1.00, (0.20, 0.50, 2)→0.80, (0.30, 0.50, 4)→1.80.

- **Coverage recommendations:** MUST test empirical VAF for representative ratios (0.25/0.50/1.0/0.0); MUST
  test Wilson center/lower/upper against the exact formula values; MUST test Wilson no-overshoot at p̂=0
  (lower=0) and p̂=1 (upper=1); MUST test `AdjustVAFForPurity` diploid het (VAF 0.4, π 0.8→1.0; VAF 0.2,
  π 0.5, n_tot 2→0.8); MUST test totalReads=0→VAF 0, altReads>totalReads→exception, negatives→exception;
  SHOULD test purity=0→exception; COULD test invariants 0≤VAF≤1 and lower≤center≤upper over a sweep.

## Deviations and assumptions

- **ASSUMPTION — Wilson `z = 1.96` for the default 95% interval.** The source states `z₀.₀₅ = 1.96`;
  used verbatim (not the more precise 1.959964) to match the cited value. Source-backed, not a free
  assumption; recorded for traceability.
- **ASSUMPTION — `AdjustVAFForPurity` fixes normal copy number = 2** and assumes the observed VAF comes
  from a heterozygous-eligible diploid-normal background (CNAqc/Tarabichi autosomal-diploid convention).
  Sex chromosomes / non-diploid normals are out of scope.

No source contradictions — GATK (empirical VAF definition), Wilson 1927 (the score interval), and CNAqc /
Tarabichi (the purity/ploidy correction) cover disjoint facets and agree.
