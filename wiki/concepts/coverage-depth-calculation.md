---
type: concept
title: "Coverage / sequencing-depth calculation (per-base depth, breadth, average)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-COVER-001-Evidence.md
  - docs/algorithms/Assembly/Coverage_Calculation.md
source_commit: 9ce49bade5c11e63eebbf8c06dd642662321d5a2
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-cover-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-COVER-001 ... Coverage (Depth) Calculation — per-base sequencing depth over a reference"
      confidence: high
      status: current
---

# Coverage / sequencing-depth calculation

Given a reference and a set of reads placed on it, **coverage (sequencing depth)** is the number
of reads that align over each reference base. This is the assembly QC read-out that answers "how
deeply is each position sequenced," and the anchor for the assembly **COVER** family. Validated
under test unit **ASSEMBLY-COVER-001**; the validation record is [[assembly-cover-001-evidence]],
and [[test-unit-registry]] tracks the unit. See [[algorithm-validation-evidence]] for the artifact
pattern.

## The three quantities (source definitions)

Traced verbatim to samtools workflows (Daniel Cook; Metagenomics Wiki) and Illumina's coverage page:

- **Per-base depth** — the number of placed reads spanning a given reference position. A position
  with no overlapping reads has depth `0`. This array is **exact**: it is a plain count, independent
  of any statistical model.
- **Average depth** — `Σ(per-base depth) / reference length` (total bases mapped / genome size).
  Equivalently the Lander-Waterman `C = LN / G` (read length × read count / genome length) for
  uniform reads.
- **Breadth of coverage** — the fraction of reference bases covered at least once:
  `(# positions with depth ≥ 1) / reference length`.

## Boundary and empty-input rules

- **Clip at the reference end** — a read extending past the reference end contributes only its
  overlapping portion to per-base depth.
- **Unmatched / empty input** — a read that does not place (below `minOverlap`) contributes `0` to
  every position; an empty read set yields an all-zero depth array of reference length, average
  depth `0` and breadth `0`.

## Worked oracle (hand-constructed, exact placement)

Reference `ACGTTGCAAT` (len 10); three distinct 5-mer reads placed unambiguously:
`ACGTT`@0 → [0,5), `TTGCA`@3 → [3,8), `GCAAT`@5 → [5,10). Depth array
`[1,1,1,2,2,2,2,2,1,1]`, Σ = 15, average = 15/10 = **1.5**, breadth = 10/10 = **1.0**.

## Lander-Waterman Poisson model (property check, not the arithmetic)

Under uniform-random read placement, per-site depth is Poisson with rate λ ≈ average depth
(Lander & Waterman 1988; Daley et al. PMC7398442). Two derived identities:

- **P(a base uncovered)** `= e^−c` (fold coverage `c`): 1× → 0.37 (~37% gaps), 5× → 0.0067.
- **Breadth** `= 1 − e^−c` (probability a base is covered ≥ once).

These are used only as a **derivation / sanity check** on the model — they are the *expected*
coverage under uniformity. The per-base depth array the unit computes is exact regardless; observed
depth can deviate from the Poisson mean when read placement is non-uniform.

## Placement is out of scope (assumption record)

The unit signature is `CalculateCoverage(reference, reads, minOverlap)`. The sources define depth
*given an alignment* but do not prescribe how an aligner places a read. The repository uses an
ungapped best-match scan requiring ≥ `minOverlap` matching characters (`FindBestAlignment`); this
placement rule is implementation-level (it decides *where* a read maps, not the counting
arithmetic). Tests therefore use exact-match reads where placement is unambiguous, isolating the
fully source-defined depth-counting rule (per-base = count of placed reads spanning the position,
clipped at the reference end). Matching is case-insensitive.
