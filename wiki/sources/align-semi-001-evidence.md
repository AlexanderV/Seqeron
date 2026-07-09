---
type: source
title: "Evidence: ALIGN-SEMI-001 (Semi-global / fitting alignment)"
tags: [validation, alignment]
doc_path: docs/Evidence/ALIGN-SEMI-001-Evidence.md
sources:
  - docs/Evidence/ALIGN-SEMI-001-Evidence.md
source_commit: c806b157357d5eccb302b3e1ea1c569f7fe48d1d
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ALIGN-SEMI-001

The validation-evidence artifact for test unit **ALIGN-SEMI-001** (Semi-Global Alignment,
fitting / query-in-reference variant). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is
summarized in [[semi-global-alignment-fitting]]. See [[test-unit-registry]] for how units
are tracked.

## What this file records

- **Online sources** (accessed 2026-02-01) — Wikipedia "Sequence alignment" and
  "Needleman–Wunsch algorithm"; Rosalind SIMS (fitting) and SMGB (semiglobal); Brudno et al.
  2003 "Glocal alignment" (doi:10.1093/bioinformatics/btg1005).
- **Algorithm spec** — semi-global as an ends-free hybrid of global + local; the fitting
  matrix initialization (first row 0, first column d·i), the NW recurrence with no zero floor,
  and traceback from the maximum of the last row.
- **Corner cases** — query embedded in / at start / at end of reference, identical sequences
  (fitting = global when lengths equal), all-mismatch (exact negative score), gap in the
  optimal path, and null input → `ArgumentNullException`.
- **Invariants (INV-1..5)** — AlignmentType = SemiGlobal; equal aligned lengths; query fully
  represented (`RemoveGaps(aligned1) == query`); reference-side substring; score = max_j F(m,j).

## Design choice (from §2.3)

The implementation deliberately selects the **query-in-reference (fitting)** member of the
semi-global family (Rosalind SIMS) rather than the overlap (OAP) or full-semiglobal (SMGB)
variants — a well-defined single mode, not a departure from the spec. No contradictions;
deviations are the fitting-variant selection plus standard .NET null-argument contract.
