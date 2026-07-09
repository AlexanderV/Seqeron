---
type: source
title: "Evidence: ALIGN-GLOBAL-001 (Needleman–Wunsch global alignment)"
tags: [validation, alignment]
doc_path: docs/Evidence/ALIGN-GLOBAL-001-Evidence.md
sources:
  - docs/Evidence/ALIGN-GLOBAL-001-Evidence.md
source_commit: 46d4efa2e08a672c942aa455eeb8b724705081e3
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ALIGN-GLOBAL-001

The validation-evidence artifact for test unit **ALIGN-GLOBAL-001** (Global Alignment,
Needleman–Wunsch). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; see that page for what every
such file contains and [[test-unit-registry]] for how units are tracked. The algorithm
itself is summarized in [[global-alignment-needleman-wunsch]].

## What this file records

- **Online sources** — Wikipedia "Needleman–Wunsch algorithm" (accessed 2026-03-06) and
  "Sequence alignment" (2026-02-01), from which the algorithm spec is traced.
- **Test dataset** — the canonical Wikipedia worked example: `GCATGCG` vs `GATTACA`,
  match +1 / mismatch −1 / gap −1, optimal score **0**, one optimal alignment
  `GCATG-CG` / `G-ATTACA` (4 matches, 2 mismatches, 2 gaps). Includes the
  border-initialized DP matrix (F(i,0) = −i, F(0,j) = −j).
- **Deviations and assumptions** — *None.* The implementation follows the standard
  linear-gap-penalty NW pseudocode exactly.

## Implementation notes (from the "Deviations" section)

- `ScoringMatrix.GapExtend` acts as the single linear gap penalty *d*; `ScoringMatrix.GapOpen`
  is **not used** by `GlobalAlign` — affine gaps are a documented NW *extension*, not part
  of the basic model.
- When several optimal alignments exist, one is returned deterministically — explicitly
  allowed by the source.
- Empty-input / null-argument handling are API-contract behaviours, not part of the NW spec.
