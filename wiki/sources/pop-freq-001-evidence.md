---
type: source
title: "Evidence: POP-FREQ-001 (Allele & genotype frequencies — MAF, MAF filtering)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-FREQ-001-Evidence.md
sources:
  - docs/Evidence/POP-FREQ-001-Evidence.md
source_commit: fec2c72b4f77c252586394fe43424909b13d98d6
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-FREQ-001

The validation-evidence artifact for test unit **POP-FREQ-001** — the population-genetics
**allele/genotype frequency** primitive: allele frequencies from genotype counts, **minor allele
frequency (MAF)** from a genotype vector, and **MAF filtering**. It is one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the formulae,
worked oracle, invariants, and corner cases are synthesized in the dedicated concept
[[allele-genotype-frequencies]]. This is a foundational `POP-*` (population-genetics) unit — the
per-population allele frequencies it produces are the fixed **F** input to
[[ancestry-estimation-admixture]] (POP-ANCESTRY-001) and the per-site `p_i` behind the
heterozygosity term of [[genetic-diversity-statistics]] (POP-DIV-001). See [[test-unit-registry]]
for how units are tracked.

## What this file records

- **Online sources (all "exact match", no deviations):**
  - **Wikipedia "Allele frequency" / "Genotype frequency"** — biallelic allele frequency from
    diploid genotype counts: `p = f(AA) + ½f(AB)`, `q = f(BB) + ½f(AB)`, **invariant `p + q = 1`**.
    Allele-counting form: `total = 2·(n_AA + n_AB + n_BB)`, `major = 2·n_AA + n_AB`,
    `minor = 2·n_BB + n_AB`. Worked oracle (four-o'clock plants) 49 AA / 42 Aa / 9 aa →
    `freq(a) = 60/200 = 0.30`, `freq(A) = 140/200 = 0.70`.
  - **Wikipedia "Minor allele frequency"** — MAF = frequency of the second-most-common allele;
    `MAF = min(alt_freq, 1 − alt_freq)`, **invariant `0 ≤ MAF ≤ 0.5`**. HapMap targeted SNPs with
    MAF ≥ 0.05; the **common (MAF > 0.05) vs rare (MAF < 0.05)** boundary, with rare variants
    enriched in coding regions.
  - **VCF/PLINK genotype dosage encoding** — 0 = homozygous ref, 1 = het, 2 = homozygous alt; alt
    count per individual = genotype value; `alt_freq = Σ genotypes / (2·n)`.
  - **Gillespie (2004) *Population Genetics: A Concise Guide*** (ISBN 978-0-8018-8008-7) and NDSU
    PopGen (McClean 1998) — textbook backing for the same formulae.

- **Datasets (documented oracles):** allele-frequency cases AF-1…AF-6 (Wikipedia 49/42/9 →
  p=0.70/q=0.30; all-hom-major → 1.0/0.0; all-hom-minor → 0.0/1.0; all-het → 0.5/0.5; HWE 25/50/25
  → 0.5/0.5; zero samples → 0,0). MAF cases MAF-1…MAF-6 (alt_freq 0.3 → MAF 0.3; alt_freq 0.7 →
  MAF 0.3; monomorphic ref/alt → 0; 50/50 → 0.5; empty → 0). Filter cases FLT-1…FLT-6 (thresholds
  MAF < 0.01, MAF < 0.05, MAF > 0.4; empty / all-filtered / none-filtered).

## Deviations and assumptions

**None** — every formula is an exact match to its source (biallelic allele frequency, `p+q=1`,
`MAF = min(f, 1−f)`, VCF dosage counting). Notable API/contract points: `CalculateAlleleFrequencies`
**validates all genotype counts ≥ 0 and throws `ArgumentOutOfRangeException` on a negative count**
(counts are non-negative by definition — validation replaces undefined behaviour); zero samples →
`(0, 0)`; empty genotype vector → MAF = 0; empty filter input → empty result. `FilterByMAF` applies
an inclusive `[minMAF, maxMAF]` band. Scope is biallelic counting/normalization only — no
Hardy–Weinberg test, no multiallelic handling, no phasing/imputation. No source contradictions;
Open Questions: none.
