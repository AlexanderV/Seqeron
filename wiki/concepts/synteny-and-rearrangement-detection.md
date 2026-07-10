---
type: concept
title: "Synteny blocks and rearrangement detection"
tags: [chromosome, comparative-genomics, algorithm]
sources:
  - docs/Evidence/CHROM-SYNT-001-Evidence.md
  - docs/algorithms/Chromosome_Analysis/Synteny_Analysis.md
  - docs/algorithms/Comparative_Genomics/Synteny_Block_Detection.md
  - docs/Validation/reports/CHROM-SYNT-001.md
  - docs/Validation/reports/COMPGEN-SYNTENY-001.md
source_commit: 3d86b2b7c044235f2082bf78748c355fefbb6176
created: 2026-07-09
updated: 2026-07-10
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
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-synteny-001-evidence
      evidence: "Test Unit ID: COMPGEN-SYNTENY-001 ... Algorithm: Synteny / Collinearity Block Detection (MCScanX collinearity model) — the comparative-genomics unit reusing this anchor"
      confidence: high
      status: current
---

# Synteny blocks and rearrangement detection

Detecting **synteny** — "conservation of blocks of order within two sets of chromosomes being
compared" (Wikipedia) — and the **chromosomal rearrangements** that break it. This is the fourth
ingested unit of the **Chromosome-analysis** family, a sibling of [[aneuploidy-detection]],
[[centromere-analysis]], and [[karyotype-analysis]]. It is deliberately positioned as the **shared
synteny anchor**: the comparative-genomics `COMPGEN-SYNTENY-001` unit covers the same syntenic-block
concept between whole genomes and should link here rather than re-deriving it. Its siblings in the
Comparative-genomics family are [[average-nucleotide-identity]] — the genome-similarity metric
(how nucleotide-identical two genomes are) to this page's gene-order-conservation view — and
[[conserved-gene-clusters-common-intervals]], the order-free counterpart that asks only whether
a gene *set* is contiguous in every genome (a common interval) rather than requiring a collinear
ordered block. The end-to-end pipeline [[genome-comparison-core-dispensable]] reuses this unit's
syntenic blocks to report an `OverallSynteny` fraction alongside its core/dispensable gene
partition. Validated under two test units: **CHROM-SYNT-001** at chromosome scale (pre-implementation evidence
[[chrom-synt-001-evidence]]; independent two-stage re-validation verdict [[chrom-synt-001-report]] —
Stage A PASS-WITH-NOTES / Stage B PASS / CLEAN, 19 tests passing, zero code change) and the
comparative-genomics **COMPGEN-SYNTENY-001** at whole-genome
scale, which reuses this anchor and supplies the concrete MCScanX DP scoring parameters (validation
record [[compgen-synteny-001-evidence]]; see the *MCScanX collinearity DP model* section below).
[[test-unit-registry]] tracks the units and [[algorithm-validation-evidence]] describes the artifact
pattern.

A **syntenic block** is "a region of chromosomes between genomes that shares a common order of
homologous genes derived from a common ancestor" (Wikipedia). The unit exposes two algorithms:
`FindSyntenyBlocks` builds blocks from ortholog pairs, and `DetectRearrangements` classifies the
structural events between adjacent blocks.

> **Note — two formulations of rearrangement detection.** This page classifies rearrangements from
> **adjacent synteny-block coordinate signals** (chromosome / strand / gap). The comparative-genomics
> [[genome-rearrangement-breakpoint-distance]] unit (COMPGEN-REARR-001) solves the same problem from a
> **signed gene-order permutation** — counting breakpoints `b(α)`, the reversal-distance lower bound
> `d≥b/2`, and classifying Inversion vs Transposition (Hannenhalli–Pevzner / Bafna–Pevzner). The two
> are complementary `alternative_to` formulations.

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

## MCScanX collinearity DP model (COMPGEN-SYNTENY-001)

The comparative-genomics unit makes the MCScanX chaining scheme behind `FindSyntenyBlocks`
explicit. Given an input **ortholog/anchor map** (anchor *generation* is delegated to
COMPGEN-ORTHO-001), a dynamic program chains adjacent anchors:

```
Score(v) = max( MatchScore(v),
                max_u [ Score(u) + MatchScore(v) + GapPenalty × NumberofGaps(u,v) ] )
MatchScore   = 50   per anchored gene pair
GapPenalty   = -1
NumberofGaps = max intervening genes between anchors u and v; must be < 25 (MAX_GAPS)
```

A non-overlapping chain is reported iff its score is **over 250** — equivalently **≥ 5 collinear
anchor pairs** (`5 × 50 = 250`), the MCScanX default minimum block. Because matches are sorted in
**both transcriptional directions**, forward and inverted (reverse) blocks are both detected
(`IsInverted`), matching the strand `+`/`-` distinction above. Anchors come from BLASTP at E ≤ 10⁻⁵
with near-duplicate collapsing (< 5 genes apart → representative pair, smallest E-value). Worked
oracles: 5 adjacent forward anchors → one forward block (score 250); reversed order → one inverted
block; 4 anchors (score 200) → no block; a ≥ 25-gene gap breaks the chain. Two source-backed
assumptions (report rule ≥ 250 **and** ≥ 5 pairs; anchors supplied as an `orthologMap`) are detailed
in [[compgen-synteny-001-evidence]]. Independent two-stage re-validation verdict
[[compgen-synteny-001-report]] — **Stage A PASS-WITH-NOTES / Stage B PASS / CLEAN**, full suite 6504
passing, no code change; the `score = n×50 − Σgaps` closed form confirmed as the cited DP recurrence's
single-chain form, and tests strengthened (+3) to lock the null-arg contract and the
direction-consistency branch. Its Stage-A notes (MAX_GAPS paper=25 vs current tool=20; "over 250"
resolved to ≥ 250 per the paper's own equivalence) are documentation-only.

## Reference tools

The definitions trace to **MCScanX** (Wang et al. 2012, synteny/collinearity detection — source of the
coordinate-only `SequenceIdentity = NaN` behaviour), **SyRI** (Goel et al. 2019, synteny and
rearrangement identification), and **MUMmer** (whole-genome alignment). No deviations from the sources
are recorded; the algorithm behaviour follows the encyclopedic synteny/rearrangement definitions.
