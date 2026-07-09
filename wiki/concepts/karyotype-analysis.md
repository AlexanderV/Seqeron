---
type: concept
title: "Karyotype analysis (chromosome-count karyotyping + ploidy detection)"
tags: [chromosome, algorithm]
sources:
  - docs/Evidence/CHROM-KARYO-001-Evidence.md
source_commit: 64c8b3182c8d371780369749eecca402938cbc44
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-karyo-001-evidence
      evidence: "Test Unit ID: CHROM-KARYO-001 ... Area: Chromosome Analysis ... AnalyzeKaryotype / DetectPloidy"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:aneuploidy-detection
      source: chrom-karyo-001-evidence
      evidence: "AnalyzeKaryotype labels aneuploidy using standard cytogenetic nomenclature based on absolute copy count (nullisomy..pentasomy), the same terminology CHROM-ANEU-001 assigns from read-depth copy number"
      confidence: high
      status: current
---

# Karyotype analysis (chromosome-count karyotyping + ploidy detection)

Building a **karyotype** — "the general appearance of the complete set of chromosomes in a cell"
(Wikipedia) — and detecting the **ploidy level** of a genome. This is the third ingested unit of the
**Chromosome-analysis** family, a sibling of [[aneuploidy-detection]] (the copy-number/read-depth anchor)
and [[centromere-analysis]] (the centromere / alpha-satellite anchor). Validated under test unit
**CHROM-KARYO-001**; the validation record is [[chrom-karyo-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

A normal human diploid karyotype is 46 chromosomes (22 autosomal pairs + 2 sex chromosomes), written
`46,XX` (female) or `46,XY` (male). The unit exposes two independent algorithms: `AnalyzeKaryotype`
karyotypes a set of chromosome **descriptors**, and `DetectPloidy` estimates ploidy from read depth.

## AnalyzeKaryotype (from chromosome descriptors)

Given a set of chromosomes (names + lengths), the algorithm:

1. Separates **sex chromosomes** from **autosomes**.
2. Groups autosomes by base chromosome name (stripping copy suffixes).
3. Counts copies of each chromosome group.
4. Compares each count against the expected ploidy level.
5. Labels aneuploidy with **standard cytogenetic nomenclature by absolute copy count**: nullisomy (0),
   monosomy (1), disomy (2), trisomy (3), tetrasomy (4), pentasomy (5) — the same ladder
   [[aneuploidy-detection]] assigns, but keyed on descriptor counts rather than depth log-ratios.

Invariants:

```
TotalChromosomes     = AutosomeCount + SexChromosomeCount
TotalGenomeSize      = Σ(chromosome lengths)
MeanChromosomeLength = TotalGenomeSize / TotalChromosomes
HasAneuploidy        = true  IFF  any chromosome group has count ≠ expectedPloidy
```

## DetectPloidy (from read depth)

Ploidy is estimated from normalized read-depth values:

1. Empty input → `(ploidy = 2, confidence = 0)` — default diploid, zero confidence.
2. Compute the **true median** of the sorted depths (average of the two middle values for even counts).
3. `ratio = medianDepth / expectedDiploidDepth`.
4. `ploidy = round(ratio × 2)`, then **clamp to `[1, 8]`**.
5. `confidence = 1.0 − |ratio × 2 − ploidy| × 2`, ∈ [0, 1].

Anchor points: diploid ratio ≈ 1.0 → ploidy 2; tetraploid ratio ≈ 2.0 → ploidy 4; haploid ratio ≈ 0.5
→ ploidy 1. Contrast [[aneuploidy-detection]]'s CN clamp `[0, 10]` and `round(2^logRatio × 2)`: here the
ratio is used directly (not through log2) and the clamp reflects a whole-genome ploidy range `[1, 8]`.

## Documented oracles

Normal diploid human (22 autosome pairs + XX/XY) → 46 chromosomes, no aneuploidy; Trisomy 21 (Down),
Monosomy X (Turner, 45,X), plus disomy/tetrasomy/pentasomy scenarios in non-diploid contexts; and the
diploid/tetraploid/haploid depth-ratio cases for `DetectPloidy`.

## Design decisions and edge cases

- **Empty chromosome input** → empty karyotype, no aneuploidy (graceful degradation).
- **Empty depth input** → `(ploidy = 2, confidence = 0)`: diploid is the common default; zero
  confidence signals "no data".
- **Ploidy clamped to `[1, 8]`.** Higher ploidy exists in nature (polytene chromosomes up to 1024-ploid
  per Wikipedia) but is outside typical analysis scope.
- **Nullisomy (0 copies) is unreachable** via `GroupBy`: the architecture detects only chromosomes
  present in the input, so an absent chromosome can never form a group — the term is mapped for
  nomenclature completeness only. A [[research-grade-limitations|research-grade]] artifact of the
  descriptor-driven design.
- **Disomy (2 copies) is aneuploidy only in non-diploid contexts.** In diploid organisms 2 copies is
  normal and never triggers; in a polyploid context 2 copies is correctly labeled Disomy per ISCN.

The artifact records **no deviations or assumptions** — the implementation follows the standard
cytogenetic definitions from the Wikipedia Karyotype / Ploidy / Aneuploidy sources.
