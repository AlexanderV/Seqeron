---
type: source
title: "Evidence: ASSEMBLY-COVER-001 (Coverage / Depth Calculation)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-COVER-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-COVER-001-Evidence.md
source_commit: 23dd89a5788feaa17e8af26acb9b4e605f6baab7
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-COVER-001

The validation-evidence artifact for test unit **ASSEMBLY-COVER-001** (per-base sequencing depth /
coverage over a reference). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm's depth/breadth/average
definitions and boundary rules are summarized in [[coverage-depth-calculation]], the anchor for the
assembly COVER family. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13, quoted verbatim):
  - **Illumina — "Sequencing Coverage for NGS Experiments"** (rank 2 vendor spec) — coverage =
    "the average number of reads that align to, or 'cover,' known reference bases"; the
    Lander-Waterman average-coverage formula **C = LN / G** (read length × read count / haploid
    genome length).
  - **Daniel E. Cook — "Calculate Depth and Breadth of Coverage From a bam File"** (rank 3 workflow)
    — depth and coverage are "the same thing, the average number of reads aligned to an individual
    base," computed as **Sum of Depths / genome size**; per-base depth via `samtools depth`; breadth
    = **Bases Mapped / genome size** (proportion of bases with ≥1 aligned read).
  - **Metagenomics Wiki — "SAMtools: get breadth of coverage"** (rank 3 tooling) — per-base depth =
    "the number of reads mapping to a specific reference position"; breadth = "fraction of reference
    bases covered by at least one read" (worked 32,876/45,678 = 0.719×); average depth = total bases
    mapped / reference length.
  - **Daley et al. 2020, PMC7398442** (rank 1 peer-reviewed) — the Lander-Waterman assumption: reads
    generated uniformly at random ⇒ per-site read count ≈ Poisson with rate λ = average depth; a site
    has "sufficient coverage" if covered by ≥ r reads.
  - **Lander & Waterman 1988, Genomics 3:231-239** (rank 1 primary / rank 4 retrieved restatement;
    full text not fetched) — gap probability **P = e^−m** (fold coverage m): 1× → 0.37, 5× → 0.0067;
    breadth complement **1 − e^−c** = probability a base is covered ≥ once.
- **Datasets** — (1) hand-constructed exact-placement oracle: reference `ACGTTGCAAT` (len 10), reads
  `ACGTT`@0 / `TTGCA`@3 / `GCAAT`@5 (distinct 5-mers, unambiguous) ⇒ depth `[1,1,1,2,2,2,2,2,1,1]`,
  Σ=15, average 15/10 = 1.5, breadth 10/10 = 1.0. (2) Lander-Waterman Poisson sanity: c=1× ⇒
  P(uncovered) e^−1 ≈ 0.3679, breadth ≈ 0.6321; c=5× ⇒ e^−5 ≈ 0.006738 — used as a property/derivation
  check only, not asserted against the exact per-base depth array.
- **Corner cases / failure modes** — position with zero overlapping reads ⇒ depth 0 (not counted in
  breadth); read extending past the reference end ⇒ only the overlapping portion contributes (boundary
  clipping); uncovered genome ⇒ breadth 0 and average depth 0; highly variable coverage — the
  uniform-Poisson mean≈C holds only under uniform random placement (a modeling caveat; the per-base
  depth array itself is exact regardless).
- **Recommended coverage** — MUST: per-base depth = count of reads spanning each position (unambiguous
  exact-match placement); read past reference end contributes only its overlap; unmatched read (below
  `minOverlap`) contributes 0 everywhere; empty reads ⇒ all-zero depth array of reference length.
  SHOULD: average depth = Σ(depth)/length = 1.5 on the worked set; breadth = (#depth≥1)/length = 1.0.
  COULD: case-insensitive matching (lowercase read vs uppercase reference still maps).

## Assumptions (from the artifact)

One assumption record: **mapping model for read placement is out of scope.** The unit signature is
`CalculateCoverage(reference, reads, minOverlap)`; the sources define depth given an *alignment* but
do not prescribe how reads are placed. The repository uses an ungapped best-match scan requiring ≥
`minOverlap` matching characters (`FindBestAlignment`) — this decides *where* a read maps, not the
depth-counting arithmetic. The counting rule itself (per-base = number of placed reads spanning the
position; clip at reference end) is fully source-defined, so tests use exact-match reads where
placement is unambiguous to isolate it.

No contradictions among the sources — Illumina, Daniel Cook, the Metagenomics Wiki, and the
Lander-Waterman papers give the same depth/average/breadth definitions (`C = LN/G`, Σdepth/G,
covered/G, breadth = 1 − e^−c). The only assumption is the out-of-scope placement model above.
