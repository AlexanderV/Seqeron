---
type: concept
title: "Runs of homozygosity & inbreeding coefficient (ROH, F_ROH)"
tags: [population-genetics, algorithm]
mcp_tools:
  - inbreeding_from_roh
  - runs_of_homozygosity
sources:
  - docs/Evidence/POP-ROH-001-Evidence.md
  - docs/algorithms/Population_Genetics/Runs_Of_Homozygosity.md
source_commit: cd2c6b2b838648c1db0897589236431a27560ecb
created: 2026-07-10
updated: 2026-07-16
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pop-roh-001-evidence
      evidence: "Test Unit ID: POP-ROH-001 ... Algorithm: Runs of Homozygosity (ROH) detection and genomic inbreeding coefficient F_ROH (methods FindROH / F_ROH)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:ancestry-estimation-admixture
      source: pop-roh-001-evidence
      evidence: "POP-ROH-001 is a population-genetics POP-* unit in the family anchored by POP-ANCESTRY-001; both operate on individual SNP genotypes but measure different quantities (per-individual homozygous-segment burden / inbreeding vs ancestry proportions)."
      confidence: high
      status: current
---

# Runs of homozygosity & inbreeding coefficient (ROH, F_ROH)

Detect the **long, uninterrupted homozygous stretches** in one individual's genome (**runs of
homozygosity**, ROH) and summarize their total burden as the **genomic inbreeding coefficient
F_ROH**. Long ROHs arise when both haplotypes descend from a recent common ancestor — the genomic
footprint of consanguinity — so their total length estimates inbreeding directly from SNP data,
without a pedigree. This is a population-genetics `POP-*` unit (**POP-ROH-001**) in the family
anchored by [[ancestry-estimation-admixture]].

It is genuinely distinct from its POP siblings: ROH is a **per-individual segment-detection** problem
over ordered genotypes along a chromosome — not per-locus counting
([[allele-genotype-frequencies]]), a within-sample diversity summary
([[genetic-diversity-statistics]]), a between-population differentiation scalar
([[population-differentiation-fst]]), a two-locus association ([[linkage-disequilibrium]]), or a
single-locus goodness-of-fit test ([[hardy-weinberg-equilibrium-test]]). Validated under test unit
**POP-ROH-001**; the literature-traced record is [[pop-roh-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

## The consecutive-runs scan (`FindROH`)

The implementation follows the **window-free consecutive-runs** method (Marras et al. 2015,
detectRUNS `consecutiveRUNS.run`) rather than PLINK's sliding window: it walks the position-sorted
SNPs one at a time and grows a run while the genotypes stay homozygous, terminating the run when a
threshold is crossed. Genotypes use the 0/1/2 allele-dosage encoding (0 = homozygous ref,
**1 = heterozygous = the "opposite" genotype**, 2 = homozygous alt), matching `CalculateMAF` /
`CalculateLD` and PLINK `--recodeA`.

A run **terminates** (and the accumulated stretch is emitted if it qualifies) when any of:

- **opposite genotypes exceed `maxOppRun`** — a few heterozygous calls are tolerated inside a run to
  absorb genotyping error; only crossing the tolerance breaks it;
- **inter-SNP gap exceeds `maxGap`** — a physical gap larger than `maxGap` breaks the run **even when
  every genotype is homozygous** (no SNP density to support the call across the gap).

A completed stretch is **retained only if it passes both** `minSNP` (minimum homozygous SNP count)
**and** `minLengthBps` (minimum physical length in bp). This dual threshold mirrors PLINK's
`--homozyg-snp` (default 100) **and** `--homozyg-kb` (default 1000 kb) — satisfying one alone is
insufficient. Each reported run carries its `Start`, `End`, and `SnpCount`.

## Inbreeding coefficient F_ROH

Once ROHs are called, the genomic inbreeding coefficient is (McQuillan et al. 2008, verbatim):

```
F_ROH = ΣL_roh / L_auto
```

`ΣL_roh` = the total length of the individual's ROHs above the chosen minimum length; `L_auto` = the
length of the SNP-covered autosomal genome (excluding centromeres) — the paper used
**2,673,768 kb (≈2,674 Mb)**. F_ROH is *the proportion of the autosomal genome sitting in long
homozygous runs*, and it tracked pedigree-derived inbreeding closely (r = 0.86 in Orkney). Worked
oracle: `ΣL_roh = 20 Mb`, `L_auto = 100 Mb` → **F_ROH = 0.20**; whole-genome coverage → 1.0.

The minimum ROH length is a **deliberate lever**: McQuillan explored ≥0.5 / ≥1.5 / ≥5 Mb — ROHs up
to ~4 Mb appear even in outbred individuals, and **1.5 Mb** best separated endogamy levels. PLINK's
sliding-window analogues (`--homozyg-window-snp 50`, `--homozyg-window-het 1`,
`--homozyg-window-threshold 0.05`, `--homozyg-density 50 kb/SNP`) parameterize the same idea a
different way.

## Invariants and edge cases

- **Each qualifying run reported once** with exact `Start ≤ End` and correct `SnpCount`.
- **One tolerated het** (`≤ maxOppRun`) keeps a **single** run; a het **beyond** tolerance **splits**
  into two runs at the correct boundaries.
- **Gap `> maxGap`** breaks an otherwise all-homozygous run.
- **Below `minSNP` or below `minLengthBps` → discarded** (both thresholds required).
- **Unsorted input is ordered internally**; **genotype `2` counts as homozygous**; **leading
  heterozygotes are skipped** (a run cannot open on an opposite genotype).
- **`F_ROH` bounded `[0, 1]`**; `0.20` and `1.0` are the traced oracle values.

## Implementation contract

Both methods live in `PopulationGeneticsAnalyzer` (`FindROH`, `CalculateInbreedingFromROH`).
A run is reported on the **inclusive `[Start, End]`** SNP positions (first/last homozygous SNP —
trailing tolerated heterozygotes are excluded); `FindROH` runs in **O(n log n)** time (dominated by
the internal position sort) and **O(n)** space over `n` SNPs. `CalculateInbreedingFromROH` is
**O(m)/O(1)** over `m` segments and computes each ROH length **half-open as `End − Start`** (so the
F_ROH `rohSegments` are `[Start, End)` intervals, distinct from the inclusive runs `FindROH` emits).

Argument validation is **eager** — the deferred `yield` iterator is wrapped so bad arguments throw
before iteration begins: `FindROH` throws `ArgumentNullException` for null genotypes and
`ArgumentOutOfRangeException` for `minSnps < 1` or negative `minLength` / `maxHeterozygotes` /
`maxGap`; `CalculateInbreedingFromROH` throws `ArgumentNullException` for null segments and returns
`0` when `genomeLength ≤ 0` (no defined denominator). ROH detection is a single linear positional
pass — not a substring search — so the repository suffix tree does not apply.

## Scope

Faithful implementation of the consecutive-runs ROH scan (Marras et al. 2015) and the
`F_ROH = ΣL_roh / L_auto` inbreeding coefficient (McQuillan et al. 2008) — no deviations. Two
documented assumptions, both API-encoding conventions: the **0/1/2 genotype encoding** (a `1` is the
opposite genotype), and **missing-genotype handling is out of scope** — `FindROH` input is
`(Position, Genotype)` with no missing sentinel, so `maxMissRun` / `--homozyg-window-missing` is not
modeled and any non-`1` genotype is treated as homozygous. It does **not** implement PLINK's exact
sliding-window scoring, missing-call bounding, LOD-based / model-based ROH callers, or ROH-island /
hotspot mapping across a cohort. No source contradictions; Open Questions: none.
