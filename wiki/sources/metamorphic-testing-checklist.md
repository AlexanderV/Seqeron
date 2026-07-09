---
type: source
title: "Checklist 02: Metamorphic Testing"
tags: [validation, testing, methodology]
doc_path: docs/checklists/02_METAMORPHIC_TESTING.md
sources:
  - docs/checklists/02_METAMORPHIC_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 02: Metamorphic Testing

The **P0** per-unit checklist for metamorphic testing — a 258-row table mapping every test unit
to the metamorphic relations (MRs) it must satisfy. Synthesized in the concept
[[metamorphic-testing]]; part of the [[validation-and-testing]] program tracked in the
[[test-unit-registry]].

## What this file records

- **Purpose:** solves the oracle problem — checks relations *between the outputs of several runs*
  under an input transformation, so no reference answer is needed.
- **Relation legend:** SUB (Subset/Superset), MON (Monotonicity), INV (Invariance), SYM
  (Symmetry), COMP (Composition), SHIFT (Positional shift).
- **Per-unit table:** all **258 units ☑ complete**, **~200+ relations** collected in
  `MetamorphicTests.cs` (started from 7 units / 18+ MRs). Representative rows: `exact ⊆
  hamming(d=0)`, prepend-flank shifts positions, `revcomp` maps a dot-plot diagonal to the
  anti-diagonal, `Fst(A,B)=Fst(B,A)`.
- **Summary:** 258 complete, ~200+ MR relations defined.

## Deviations and contradictions

None. The relation legend deliberately overlaps the [[property-based-testing]] invariant legend
(INV/SYM/MON vs I/S/M) — the same mathematical structure applied across runs rather than within
one run.
