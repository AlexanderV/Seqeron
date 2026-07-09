---
type: source
title: "Evidence: CHROM-SYNT-001 (Synteny analysis — synteny blocks + rearrangement detection)"
tags: [validation, chromosome]
doc_path: docs/Evidence/CHROM-SYNT-001-Evidence.md
sources:
  - docs/Evidence/CHROM-SYNT-001-Evidence.md
source_commit: ba8861eae3de465dbda5246943d5465da6af6389
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CHROM-SYNT-001

The validation-evidence artifact for test unit **CHROM-SYNT-001** — synteny analysis: finding syntenic
blocks from ortholog pairs (`FindSyntenyBlocks`) and classifying the chromosomal rearrangements between
them (`DetectRearrangements`). This is the fourth **Chromosome-analysis** family Evidence file (after
[[chrom-aneu-001-evidence]], [[chrom-cent-001-evidence]], and [[chrom-karyo-001-evidence]]) and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the two algorithms, their invariants, the documented oracles, and the edge cases are summarized in
[[synteny-and-rearrangement-detection]], the shared synteny anchor (the comparative-genomics
`COMPGEN-SYNTENY-001` unit will reuse that concept). See [[test-unit-registry]] for how units are
tracked.

## What this file records

- **Online sources:**
  - **Wikipedia "Synteny"** — synteny commonly means colinearity: conservation of blocks of gene order
    across two chromosome sets; these regions are **syntenic blocks**.
  - **Wikipedia "Comparative genomics"** — a syntenic block is a chromosomal region sharing a common
    order of homologous genes derived from a common ancestor.
  - **Wikipedia "Chromosomal rearrangement"** — the four event types: inversion (segment reversed),
    translocation (segment to a different chromosome), deletion (segment removed), duplication (segment
    copied).
  - **Primary literature:** MCScanX (Wang et al. 2012, *Nucleic Acids Res.* 40(7):e49), SyRI (Goel et
    al. 2019, *Genome Biology* 20:277), Liu et al. 2018 (*BMC Bioinformatics* 19(1):26, systematic
    synteny-inference evaluation). Tool references: MCScanX, SyRI, MUMmer.
- **Algorithm behaviour (from the artifact):**
  - **`FindSyntenyBlocks`** — group ortholog pairs by chromosome pair; sort by reference position;
    identify collinear runs; merge consecutive segments respecting `maxGap`; emit blocks meeting
    `minGenes`. Params `minGenes` (default 3–5) and `maxGap` (megabases). Invariants I1–I5:
    `GeneCount ≥ minGenes`, valid coordinates (`Start ≤ End`) for both species, `Strand ∈ {'+','-'}`,
    `SequenceIdentity = NaN` (not computable from coordinate-only input per MCScanX), all genes in one
    chromosome pair.
  - **`DetectRearrangements`** — sort blocks by reference chr/position; compare adjacent blocks:
    different target chr → **Translocation**, same target chr + different strand → **Inversion**, gap
    asymmetry → **Deletion**, overlapping source coords + different targets → **Duplication**.
    Invariants: `Type` is a recognized value, `Position1` always non-null, and for translocations
    `Chromosome2` differs from the source block's target chromosome.
- **Datasets (documented oracles):** collinear forward block (4 genes chr1→chrA → 1 block `+`,
  GeneCount 4, 1000–8000); inverted block (reverse order → 1 block `-`); translocation detection
  (chr1→chrA then chr1→chrB → translocation at 50000); inversion detection (chr1→chrA `+` then `-` →
  inversion, positions 50000/60000, size 10000).

## Coverage classification changes (artifact §7, 2026-03-08)

| Action | Count | Details |
|--------|-------|---------|
| ⚠ Weak → Strengthened | 8 | M1, M2, M5, M6, M9, M10, M14, M15: range/permissive assertions replaced with exact hand-calculated values |
| 🔁 Duplicate → Removed | 2 | M7 (GeneCountMatchesInput), M8 (CoordinatesSpanAllGenes): subsumed by strengthened M1 |
| ❌ Missing → Implemented | 1 | M16 (GapExceedsMaxGap_SplitsIntoSeparateBlocks): the `maxGap` parameter was untested |

## Deviations and assumptions

The artifact records **no deviations** — the implementation follows the encyclopedic synteny and
rearrangement definitions, with the MCScanX-backed `SequenceIdentity = NaN` for coordinate-only input.
The edge-case rows sourced as "Implementation" (empty input → empty) are standard API-contract
behaviour outside the algorithm spec. No contradictions among sources: the Wikipedia synteny /
comparative-genomics / rearrangement definitions and the MCScanX/SyRI tool descriptions agree.
