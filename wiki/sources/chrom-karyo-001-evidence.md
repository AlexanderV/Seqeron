---
type: source
title: "Evidence: CHROM-KARYO-001 (Karyotype analysis — karyotyping + ploidy detection)"
tags: [validation, chromosome]
doc_path: docs/Evidence/CHROM-KARYO-001-Evidence.md
sources:
  - docs/Evidence/CHROM-KARYO-001-Evidence.md
source_commit: 64c8b3182c8d371780369749eecca402938cbc44
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CHROM-KARYO-001

The validation-evidence artifact for test unit **CHROM-KARYO-001** — karyotype analysis: karyotyping a
set of chromosome descriptors (`AnalyzeKaryotype`) and detecting the genome's ploidy level from read
depth (`DetectPloidy`). This is the third **Chromosome-analysis** family Evidence file (after
[[chrom-aneu-001-evidence]] and [[chrom-cent-001-evidence]]) and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the two algorithms, their
invariants, the documented oracles, and the five design decisions are summarized in
[[karyotype-analysis]], the anchor concept. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (all Wikipedia, verified 2026-03-08):
  - **"Karyotype"** — a karyotype is "the general appearance of the complete set of chromosomes";
    normal human diploid = 46 (22 autosomal pairs + 2 sex chromosomes); notation `46,XX` / `46,XY`;
    autosomes numbered 1–22 largest-to-smallest; human chromosome groups A–G by size and centromere
    position.
  - **"Ploidy"** — number of complete chromosome sets: monoploid (1) / diploid (2) / triploid (3) /
    tetraploid (4) / …; humans 2n = 46, n = 23; euploidy (normal set count) vs aneuploidy (abnormal
    individual chromosome count); polyploidy common in plants, rare in animals.
  - **"Aneuploidy"** — abnormal count of *individual* chromosomes (not whole-set); terminology by
    **absolute copy count**: nullisomy 0 / monosomy 1 (Turner 45,X) / disomy 2 (normal diploid) /
    trisomy 3 (Down, Trisomy 21) / tetrasomy 4 / pentasomy 5; survivable autosomal trisomies 21/18/13.
- **Algorithm behaviour (from the implementation):**
  - **`AnalyzeKaryotype`** — separate sex chromosomes from autosomes; group autosomes by base name
    (strip copy suffixes); count copies per group; compare against expected ploidy; label aneuploidy
    with standard cytogenetic nomenclature by absolute copy count. Invariants:
    `TotalChromosomes = AutosomeCount + SexChromosomeCount`, `TotalGenomeSize = Σ(lengths)`,
    `MeanChromosomeLength = TotalGenomeSize / TotalChromosomes`, `HasAneuploidy` iff any group count ≠
    expected ploidy.
  - **`DetectPloidy`** — empty input → `(2, 0)`; compute the true median of sorted depths (average of
    two middle values for even counts); `ratio = medianDepth / expectedDiploidDepth`;
    `ploidy = round(ratio × 2)` clamped to `[1, 8]`; `confidence = 1.0 − |ratio × 2 − ploidy| × 2`.
    Invariants: PloidyLevel ∈ [1, 8], Confidence ∈ [0, 1]; diploid ratio ≈ 1.0, tetraploid ≈ 2.0.
- **Datasets (documented oracles):** normal diploid human (22 pairs + XX/XY → 46, no aneuploidy);
  Trisomy 21 (Down, 3× chr21 → "Trisomy"); Turner (45,X, single X → "Monosomy"); disomy in a
  tetraploid; tetrasomy / pentasomy of a chromosome in a diploid; and diploid (ratio ≈ 1.0 → 2),
  tetraploid (≈ 2.0 → 4), haploid (≈ 0.5 → 1) depth cases.

## Design decisions (from the artifact §"Design Decisions")

| ID | Decision | Rationale |
|----|----------|-----------|
| DD1 | Empty chromosome input → empty karyotype, no aneuploidy | Graceful degradation |
| DD2 | Empty depth input → `(ploidy=2, confidence=0)` | Diploid default; zero confidence = no data |
| DD3 | Ploidy clamped to `[1, 8]` | Practical limit (up to 1024-ploid exists in nature, out of scope) |
| DD4 | Nullisomy (0 copies) unreachable via `GroupBy` | Absent chromosomes form no group; term mapped for completeness |
| DD5 | Disomy (2 copies) is aneuploidy only in non-diploid contexts | 2 copies is normal in diploids; correctly labeled Disomy per ISCN in polyploid contexts |

## Deviations and assumptions

The artifact's "Deviations and Assumptions" section is **None** — the implementation follows the
standard cytogenetic definitions from the three Wikipedia sources exactly. No contradictions among the
sources; DD4/DD5 are architecture/nomenclature notes, not departures from the spec.
