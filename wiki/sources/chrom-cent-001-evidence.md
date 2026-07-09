---
type: source
title: "Evidence: CHROM-CENT-001 (Centromere analysis — arm-ratio class + alpha-satellite/HOR)"
tags: [validation, chromosome]
doc_path: docs/Evidence/CHROM-CENT-001-Evidence.md
sources:
  - docs/Evidence/CHROM-CENT-001-Evidence.md
source_commit: a973e9b9fb6201395fd56c6177eff279548669e2
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CHROM-CENT-001

The validation-evidence artifact for test unit **CHROM-CENT-001** — centromere analysis:
`ChromosomeAnalyzer.AnalyzeCentromere(...)` plus four opt-in, additive detectors layered onto it over
successive sessions. This is the second **Chromosome-analysis** family Evidence file (after
[[chrom-aneu-001-evidence]]) and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the Levan arm-ratio classification, the
alpha-satellite / CENP-B / HOR / suprachromosomal-family detection layers, the derived parameters and
their limits are synthesized in [[centromere-analysis]], the anchor for the chromosome
centromere/satellite family. See [[test-unit-registry]] for how units are tracked. Unusually for this
family, the file is a **layered, multi-session** record (base + four dated addenda), not a single
templated pass.

## What this file records

- **Online / literature sources:**
  - **Wikipedia — "Centromere", "Karyotype", "Chromosome"** (encyclopedia) — centromere links sister
    chromatids, defines the short (p, "petit") and long (q) arms; the five position classes; human
    karyotype groups A/D/G by centromere position; constitutive heterochromatin of repetitive DNA.
  - **Levan, Fredga & Sandberg (1964)**, *Hereditas* 52(2):201–220 — the arm-ratio (q/p) nomenclature
    for centromeric position (the classification-threshold source, cited via Wikipedia).
  - **Hartley & O'Neill (2019)** / **McNulty & Sullivan (2018)**, *Genes / Chromosome Res*
    (PMC6121732) — 171-bp alpha-satellite monomer; 50–70% intra-HOR monomer identity; HOR-period
    definition; 97–100% inter-HOR identity within an array; A/B-type monomer taxonomy.
  - **Masumoto et al. (1989)**, *J Cell Biol* 109(4):1963 (via PMC4843215) — the 17-bp CENP-B box, exact
    IUPAC consensus `YTTCGTTGGAARCGGGA`.
  - **Rosandić, Paar et al. (2024)** (PMC11050224) + **Alkan et al. 2007 / ColorHOR (Paar 2005)**
    (Oxford Academic) — n-mer HOR cascade, <5% inter-HOR divergence, the chr1 1866-bp 11-mer worked
    example (168–171-bp monomers, 1.8% divergence).
  - **Shepelev et al. (2009)**, *PLoS Genet* 5(9):e1000641 + **Dfam** ALR/ALRa/ALRb consensus monomers
    (CC0) + **T2T/CHM13** (CC0) — the suprachromosomal-family (SF1–SF5) taxonomy and the bundled
    CC0 reference.

- **Algorithm behaviour (from the implementation):**
  - **Base `AnalyzeCentromere`** — sliding-window k-mer frequency, GC-variability discriminator
    (centromeres have low GC variability), repeat-content estimate (k=15). `AlphaSatelliteContent` is a
    **generic tandem-repeat-density** score, *not* an alpha-satellite-specific measurement.
  - **Levan arm-ratio classification** (exact to the 1964 table): ratio ≤ 1.7 → Metacentric;
    (1.7, 3.0] → Submetacentric; (3.0, 7.0) → Subtelocentric; ≥ 7.0 → Acrocentric; p = 0 → Telocentric.
  - **`DetectAlphaSatellite` / `FindCenpBBoxes`** (added 2026-06-24, additive) — alpha-satellite-specific
    signal from tandem 171-bp periodicity (±5 bp tolerance, ≥0.50 self-similarity) + AT-richness (>0.50)
    + IUPAC CENP-B box match; **no consensus monomer string is embedded**.
  - **`DetectHigherOrderRepeat`** (added 2026-06-24, opt-in) — split the array into 171-bp monomers,
    align monomer-vs-monomer with `SequenceAligner.GlobalAlign` + `CalculateStatistics`; HOR period =
    smallest block size k with monomers k apart ≥95% identical (inter-HOR <5% divergence) across ≥90%
    of the array; period 1 = homogeneous 1-mer array (not a multi-monomer HOR).
  - **`AssignSuprachromosomalFamily`** (added 2026-06-25, opt-in) — best-match each monomer to a
    **bundled CC0 Dfam reference** (ALR/ALRa = A-type, ALRb = B-type, ≥60% identity gate), take the HOR
    period, assign SF from period + A/B composition (SF3 pentameric period%5==0, SF4 monomeric A-type,
    {SF1,SF2} dimeric A→B, SF5 irregular).

- **Datasets / oracles:** the Levan threshold table (M/m/sm/st/t/T signs at ratios 1.7/3.0/7.0/∞);
  human chromosome centromere positions (chr1 metacentric, chr2 submetacentric, chr13/21/Y acrocentric);
  the chr1 1866-bp 11-mer HOR with 1.8% inter-copy divergence.

## Assumptions / limitations (from the artifact)

- The `AlphaSatelliteContent` field of the base method is a **repeat-density heuristic**, not an
  alpha-satellite-specific measurement; the specific signal lives in the additive detectors.
- Practical human karyotype classes come from **cytogenetic (microscopic) observation**, not genomic
  coordinate ratios; Levan thresholds on sequence-derived positions may reclassify borderline
  chromosomes.
- Two derived parameters are flagged **ASSUMPTION**: the ≥60% identity gate to call alpha-satellite
  vs a reference, and the SF3 ⇔ period%5==0 pentameric-ancestry proxy.
- **Residual, data-blocked:** the bundled CC0 reference resolves SF3/SF4/SF5 and narrows dimeric arrays
  to {SF1, SF2} but does **not** separate SF1 from SF2 (both dimeric with identical A→B box pattern) nor
  tag diverged-pentamer SF3 arrays whose period isn't a multiple of 5 (e.g. DXZ1 period 12). Doing
  either needs an SF-resolved consensus monomer library (J1/J2/D1/D2/W1–W5/M1/R1/R2), published only in
  non-machine-retrievable supplements or unlicensed HMM repos — **not redistributable** (cf. the
  TIGRFAM/LM22 non-redistribution rule). Callers holding such a set can pass it to
  `AssignSuprachromosomalFamily(sequence, reference)`.

No contradictions among the sources — the encyclopedic articles, Levan nomenclature, the alphoid-DNA
literature (Hartley/McNulty/Masumoto/Rosandić/Alkan/Shepelev), and the Dfam/T2T reference agree; the
171-bp monomer, 17-bp CENP-B box, and <5% inter-HOR divergence recur consistently across them.
