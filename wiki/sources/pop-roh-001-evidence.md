---
type: source
title: "Evidence: POP-ROH-001 (Runs of homozygosity + inbreeding coefficient F_ROH)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-ROH-001-Evidence.md
sources:
  - docs/Evidence/POP-ROH-001-Evidence.md
source_commit: 0885cd52a7dfdf9158d2607127ceea1945c8dff2
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-ROH-001

The validation-evidence artifact for test unit **POP-ROH-001** — detection of **runs of
homozygosity (ROH)** (`FindROH`) and the genomic inbreeding coefficient **F_ROH**. It is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the method, formula, invariants, worked oracles, and corner cases are synthesized in the dedicated
concept [[runs-of-homozygosity-inbreeding]]. This is a population-genetics `POP-*` unit in the family
anchored by [[ancestry-estimation-admixture]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **McQuillan et al. (2008)** *Runs of Homozygosity in European Populations*, Am J Hum Genet
    83(3):359–372 (rank 1, primary). Verbatim **F_ROH = ΣL_roh / L_auto** — proportion of the
    autosomal genome in ROHs above a specified minimum length; `ΣL_roh` = total length of an
    individual's ROHs, `L_auto` = length of SNP-covered autosomal genome (excludes centromeres),
    with the numeric value **2,673,768 kb (≈2,674 Mb)** used in the paper. F_ROH correlated with
    pedigree inbreeding in Orkney (r = 0.86). Length thresholds explored: ≥0.5 / ≥1.5 / ≥5 Mb (ROHs
    up to ~4 Mb occur in outbred individuals; 1.5 Mb best distinguished endogamy levels).
  - **Chang et al. (2015)** PLINK 1.9 `--homozyg` documentation (rank 3, reference impl). The
    sliding-window default parameters: `--homozyg-snp 100`, `--homozyg-kb 1000` (both a **SNP-count
    AND a physical-length threshold**, both must pass), `--homozyg-window-snp 50`,
    `--homozyg-window-het 1`, `--homozyg-window-missing 5`, `--homozyg-window-threshold 0.05`,
    `--homozyg-gap 1000 kb`, `--homozyg-density 50 kb/SNP`.
  - **Marras et al. (2015)** consecutive-runs method / detectRUNS `consecutiveRUNS.run` (rank 3,
    reference impl of Anim Genet 46(2):110–121). The **window-free** scan the implementation follows:
    walk the genome SNP by SNP; a run continues while criteria hold and terminates when a threshold
    is violated. Parameters `maxOppRun` (max opposite/heterozygous genotypes tolerated in a run),
    `maxMissRun` (max missing), `minSNP` (min homozygous SNPs to retain), `minLengthBps` (min
    physical length), `maxGap` (max inter-SNP distance before the run breaks).

- **Corner cases / failure modes (documented):**
  - A small number of heterozygous calls (`≤ maxOppRun`) is tolerated inside a run to absorb
    genotyping error; only *exceeding* the tolerance breaks the run.
  - An inter-SNP gap larger than `maxGap` breaks a run even when every genotype is homozygous.
  - Stretches failing `minSNP` **or** `minLengthBps` are discarded (PLINK: both count and length
    thresholds must be satisfied — passing one alone is insufficient).

- **Oracles / datasets:**
  - *McQuillan F_ROH constants* — `L_auto = 2,673,768 kb`; worked example `ΣL_roh = 20 Mb`,
    `L_auto = 100 Mb` → **F_ROH = 0.20**; whole-genome ΣL_roh = L_auto → F_ROH = 1.0.
  - *PLINK default thresholds* — min 100 SNPs, min 1000 kb length, ≤1 het/window, max gap 1000 kb.

## Deviations and assumptions

**No algorithm deviations** — the consecutive-runs termination logic (Marras 2015) and the
`F_ROH = ΣL_roh / L_auto` formula (McQuillan 2008) are followed as specified. Two documented
**assumptions**, both API-encoding notes rather than invented behaviour:

1. **Genotype encoding 0/1/2** — 0 = homozygous reference, 1 = heterozygous, 2 = homozygous
   alternate (the additive allele-dosage convention shared by `CalculateMAF`/`CalculateLD` and PLINK
   `--recodeA`). A `1` is therefore the *opposite* (heterozygous) genotype counted against
   `maxOppRun`; the rule is encoding-independent (homozygous vs opposite).
2. **Missing-genotype handling out of scope** — `FindROH` input is `(Position, Genotype)` with no
   missing sentinel, so `maxMissRun` / `--homozyg-window-missing` is not modeled; any non-`1`
   genotype is treated as homozygous. Documented limitation, not invented behaviour.

Recommended coverage (MUST): single uninterrupted run reported once with exact Start/End/SnpCount;
one tolerated het (≤ maxOppRun) keeps one run; a het beyond tolerance splits into two runs at correct
boundaries; gap > maxGap breaks an all-homozygous run; runs below minSnps or minLength discarded;
**F_ROH = ΣL_roh / L_auto** exact (0.20 example; 1.0 whole-genome). SHOULD: unsorted input ordered
internally, genotype 2 counts as homozygous, leading heterozygotes skipped. COULD: invalid arguments
(null, negative thresholds, non-positive genome length). No contradictions among sources; Open
Questions: none.
