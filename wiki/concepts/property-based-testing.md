---
type: concept
title: "Property-based testing (FsCheck) — invariants over generated inputs"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/01_PROPERTY_BASED_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: property-based-testing-checklist
      evidence: "Priority P0; one of the ten complementary testing methodologies. 'Кожен геномний алгоритм має щонайменше один виражуваний інваріант' — every algorithm carries at least one expressible invariant, checked over hundreds of generated inputs via FsCheck.NUnit."
      confidence: high
      status: current
---

# Property-based testing (FsCheck)

Property-based testing generates **hundreds of random inputs** and asserts that the algorithm
satisfies a **mathematical invariant** on every one, rather than checking a handful of
hand-picked input→output pairs. The generator explores the input space; the property is the
oracle. Seqeron uses **FsCheck** wired into NUnit via `FsCheck.NUnit`, and the discipline is
that **every genomic algorithm must express at least one invariant**. This is the **P0**
(highest-priority) member of the ten-methodology [[validation-and-testing]] program, the
complement that the [[advanced-testing-checklist]] flagged as the biggest early gap. The
checklist record is [[property-based-testing-checklist]].

## The invariant taxonomy

Seven reusable invariant shapes cover the library; a property file names which apply to each unit:

- **R — Range**: the result lands in its admissible interval (`GC% ∈ [0,100]`, `Fst ∈ [0,1]`,
  edit distance `≥ 0`).
- **S — Symmetry**: `f(a,b) = f(b,a)` (alignment score, Hamming/edit distance, Fst, ANI).
- **I — Idempotence / Involution**: `f(f(x)) = x` — `revcomp(revcomp(x)) = x`,
  `complement(complement(x)) = x`, `d(x,x) = 0`.
- **M — Monotonicity**: more of X drives Y up or down (lower `minRepeats` → ≥ results; more
  GC → higher Tm; more mismatches → lower off-target score).
- **P — Preservation**: a property survives a transformation (complement preserves GC%; length
  preserved under revcomp; per-amino-acid codon fractions sum to 1).
- **RT — Round-trip**: `parse(serialize(x)) = x` — the parser/serializer isomorphism for FASTA,
  FASTQ, GFF, GenBank, EMBL, Newick, dot-bracket.
- **D — Determinism**: same input → same output (asserted almost everywhere).

## Coverage

**258 / 258 test units complete** at the 2026-03-19 checklist. ~22 files under `Properties/`
(GcContent, Sequence, EditDistance, Alignment, Kmer, Phylogenetic, PopGen, PrimerProbe, RNA,
Splicing, …), extended with new files for the Chromosome, Epigenetics, Oncology, Statistics,
Comparative, Assembly, PanGenome, Transcriptome and StructuralVariant families as those units
landed. This is the load-bearing layer of the testing program: it is what makes the
[[research-grade-limitations|research-grade]] correctness claim more than assertion.

## Where it sits among siblings

Property-based testing is the general form; two siblings sharpen it. [[algebraic-testing]]
verifies the *named laws* (identity, associativity, involution, round-trip) more formally, and
[[metamorphic-testing]] handles the case where **no oracle exists** — relating outputs of
multiple runs instead of checking one. The FsCheck invariants and the metamorphic relations
overlap deliberately (both encode symmetry, monotonicity, involution); the difference is that a
property asserts on one run while a metamorphic relation ties several runs together. Each unit
also clears the [[definition-of-done]] and is tracked in the [[test-unit-registry]].
