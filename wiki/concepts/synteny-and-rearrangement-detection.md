---
type: concept
title: "Synteny blocks and rearrangement detection"
tags: [chromosome, comparative-genomics, algorithm]
sources:
  - docs/Evidence/CHROM-SYNT-001-Evidence.md
source_commit: ba8861eae3de465dbda5246943d5465da6af6389
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-synt-001-evidence
      evidence: "Test Unit ID: CHROM-SYNT-001 ... Area: Chromosome Analysis ... FindSyntenyBlocks / DetectRearrangements"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:aneuploidy-detection
      source: chrom-synt-001-evidence
      evidence: "CHROM-SYNT-001 is a sibling Chromosome-analysis unit of CHROM-ANEU-001; both operate at whole-chromosome scale (synteny/rearrangement vs copy-number)"
      confidence: medium
      status: current
---

# Synteny blocks and rearrangement detection

Detecting **synteny** — "conservation of blocks of order within two sets of chromosomes being
compared" (Wikipedia) — and the **chromosomal rearrangements** that break it. This is the fourth
ingested unit of the **Chromosome-analysis** family, a sibling of [[aneuploidy-detection]],
[[centromere-analysis]], and [[karyotype-analysis]]. It is deliberately positioned as the **shared
synteny anchor**: the comparative-genomics `COMPGEN-SYNTENY-001` unit covers the same syntenic-block
concept between whole genomes and should link here rather than re-deriving it. Validated under test
unit **CHROM-SYNT-001**; the validation record is [[chrom-synt-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

A **syntenic block** is "a region of chromosomes between genomes that shares a common order of
homologous genes derived from a common ancestor" (Wikipedia). The unit exposes two algorithms:
`FindSyntenyBlocks` builds blocks from ortholog pairs, and `DetectRearrangements` classifies the
structural events between adjacent blocks.

## FindSyntenyBlocks (blocks from ortholog pairs)

Input: a list of **ortholog pairs** (each with a position in both genomes) plus two parameters —
`minGenes` (minimum genes to form a block, typical default 3–5) and `maxGap` (maximum allowed gap
between consecutive genes, in megabases). The algorithm:

1. Groups ortholog pairs by **chromosome pair** (reference chr → target chr).
2. Sorts by position in the reference genome.
3. Identifies **collinear runs** — consecutive genes keeping the same relative order.
4. Merges consecutive collinear segments while respecting the `maxGap` constraint.
5. Emits only blocks meeting the `minGenes` threshold.

Each block carries coordinates in both species, a **strand** (`+` forward/same order, `-`
inverted/opposite order), a gene count, and a sequence-identity field.

Invariants:

```
GeneCount        >= minGenes            (I1)
Start <= End     for both species       (I2, valid coordinates)
Strand           in { '+', '-' }        (I3)
SequenceIdentity = NaN                  (I4, not computable from coordinate-only input, per MCScanX)
all genes        in one chromosome pair (I5)
```

## DetectRearrangements (events between blocks)

Input: a list of synteny blocks. The algorithm sorts blocks by reference chromosome and position,
then compares **adjacent** blocks to classify the rearrangement type (Wikipedia
"Chromosomal rearrangement"):

| Adjacent-block signal | Event |
|-----------------------|-------|
| Different target chromosome | **Translocation** |
| Same target chromosome, different strand | **Inversion** |
| Same chr pair / strand, gap asymmetry | **Deletion** |
| Overlapping source coords, different targets | **Duplication** |

Invariants:

```
Type        in { "Inversion", "Translocation", "Deletion", "Duplication" }   (I1)
Position1   always set (non-null)                                             (I2)
Chromosome2 differs from source block's target chromosome (for translocations) (I3)
```

## Documented oracles

- **Collinear forward block** — 4 genes chr1→chrA in forward order → 1 block, strand `+`, GeneCount 4,
  spans 1000–8000 in both species.
- **Inverted block** — 4 genes in reverse order in the target → 1 block, strand `-`, GeneCount 4.
- **Translocation** — Block1 chr1→chrA, Block2 chr1→chrB → translocation (chrA→chrB) at position 50000.
- **Inversion** — Block1 chr1→chrA `+`, Block2 chr1→chrA `-` → inversion, positions 50000/60000,
  size 10000.

## Edge cases

- **Empty input** → empty result (both algorithms).
- **Fewer genes than `minGenes`** (single gene, or two genes with `minGenes=3`) → empty.
- **Gap exceeds `maxGap`** → the run breaks into separate blocks (the M16 coverage gap this unit
  closed — `maxGap` was previously untested).
- **Multiple chromosome pairs** → separate blocks per pair.
- **Single block** into `DetectRearrangements` → empty (no adjacent pairs to compare).
- **All same chromosome, same strand** → empty (no rearrangement).

## Reference tools

The definitions trace to **MCScanX** (Wang et al. 2012, synteny/collinearity detection — source of the
coordinate-only `SequenceIdentity = NaN` behaviour), **SyRI** (Goel et al. 2019, synteny and
rearrangement identification), and **MUMmer** (whole-genome alignment). No deviations from the sources
are recorded; the algorithm behaviour follows the encyclopedic synteny/rearrangement definitions.
