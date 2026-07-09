---
type: concept
title: "Metamorphic testing — relations across runs when no oracle exists"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/02_METAMORPHIC_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: metamorphic-testing-checklist
      evidence: "Priority P0. 'Metamorphic testing розв'язує проблему оракула' — solves the oracle problem by checking metamorphic relations that link the outputs of several runs under an input transformation, instead of checking a single expected result."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:property-based-testing
      source: metamorphic-testing-checklist
      evidence: "Both encode symmetry/monotonicity/involution invariants; the metamorphic relation ties several runs together (e.g. exact ⊆ approximate, prepend-flank shifts positions) whereas a property asserts on one run. The checklist legend (INV/SYM/MON/SUB/COMP/SHIFT) overlaps the property taxonomy (R/S/I/M/P/RT/D)."
      confidence: high
      status: current
---

# Metamorphic testing

Metamorphic testing attacks the **oracle problem** — the situation where, for an arbitrary
input, you have no reference answer to compare against (what *is* the "correct" MSA, the correct
folded structure, the correct phylogeny?). Instead of asserting a specific output, it asserts a
**metamorphic relation (MR)**: a property that must hold *between the outputs of two or more
runs* whose inputs are related by a known transformation. If `revcomp` maps a diagonal dot-plot
to an anti-diagonal, or prepending a flank shifts every reported position by the flank length,
the relation is checkable without ever knowing the absolute answer. This is a **P0** member of
the [[validation-and-testing]] program and the natural partner to [[property-based-testing]];
the checklist record is [[metamorphic-testing-checklist]].

## The relation taxonomy

- **SUB — Subset/Superset**: widening a parameter yields a superset (`exact ⊆ hamming(d=0)`;
  lower score threshold → ≥ matches; `N ⊇` all four bases in IUPAC degeneracy).
- **MON — Monotonicity**: more of X → more/less Y across runs (add GC → Tm up; more mismatches
  → lower off-target score; more segregating sites → higher θ).
- **INV — Invariance**: a transformation of the input leaves the output unchanged (shuffle
  preserves GC%; permuting rows preserves a consensus; adding a distant flank doesn't change a
  local alignment).
- **SYM — Symmetry**: `f(a,b) = f(b,a)` across the swapped run (Fst, RF distance, ANI).
- **COMP — Composition**: outputs of two operations compose or nest
  (`parse(write(x)) = x`; local alignment of identical sequences = global; gene ⊃ ORF).
- **SHIFT — Positional shift**: prepending a flank shifts every reported coordinate by the
  flank length (motifs, ORFs, repeats, restriction sites, breakpoints).

## Coverage

**258 / 258 units** at the 2026-03-19 checklist, **~200+ relations** collected in
`MetamorphicTests.cs`. The relations are especially valuable for the search/alignment/scoring/
statistics families where a golden answer is unavailable but the transformation behaviour is
mathematically pinned. Together with [[property-based-testing]] and [[algebraic-testing]] this
forms the invariant-driven core of the program; [[snapshot-testing]] covers the complementary
case where the exact output *is* fixed and any diff must be reviewed.
