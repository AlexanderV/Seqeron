---
type: source
title: "Evidence: RNA-PSEUDOKNOT-001 (Pseudoknot detection — crossing base pairs in a given structure)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md
sources:
  - docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md
source_commit: ae0dfc54b6b719a8fa68c2f120f3f4e3235cd02e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-PSEUDOKNOT-001

The validation-evidence artifact for test unit **RNA-PSEUDOKNOT-001** — **pseudoknot detection**:
given a **set of base pairs** (a structure that is already known), identify the **crossing** pairs
that make it pseudoknotted. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept is
[[rna-pseudoknot-detection]]. [[test-unit-registry]] tracks the unit.

This is the **detection / analysis** counterpart of the crossing-helix family, distinct from the
energy-driven [[rna-pseudoknot-prediction|pknotsRG predictor]]: **detection takes a base-pair set as
input and reports crossings** (a pure O(n²) combinatorial scan, no sequence and no energy model),
whereas prediction takes a **sequence** and *folds* it (O(n⁴), Turner NN energy). Detection is exactly
the primitive the predictor's validity invariant leans on (`DetectPseudoknots` finds ≥1 genuine
crossing when a knot is returned).

## What this file records

- **Algorithm:** `RnaSecondaryStructure.DetectPseudoknots(basePairs)` — scan a set of base pairs for
  every pair-of-pairs that **cross**, reporting each crossing as one pseudoknot.

- **Crossing condition (the definitional core, sourced verbatim):** two pairs cross when, renaming to
  this unit's symbols, **`i < k < j < l`** — one pair's opener lies inside the other pair while its
  closer lies outside (Antczak et al. 2018: "for any pair (*i, i'*) there exists another one (*j, j'*)
  such that *i < j < i' < j'*"). The two exhaustive negatives are **nested** (`i < k < l < j`, one pair
  fully inside the other) and **disjoint** (`j < k`, side-by-side ranges) — neither is a pseudoknot.

- **Authoritative sources:**
  - **Antczak et al. (2018)** *Bioinformatics* 34(8):1304 (rank 1) — the crossing/conflict definition
    `i < k < j < l`, pseudoknot **order** (minimum number of base-pair-set decompositions yielding a
    nested structure), and the dot-bracket-letter (DBL) notation by region order
    (order 0 `()`, 1 `[]`, 2 `{}`, 3 `<>`, 4–8 letters); H-type knot = `([)]`.
  - **Smit, Rother, Heringa & Knight (2008)** *RNA* 14(3):410 (rank 1) — pseudoknot presence *requires*
    crossing pairs; obtaining a nested structure means disregarding some pairs; origin of the
    order-assignment / pseudoknot-removal family biotite implements.
  - **biotite.structure.pseudoknots** (rank 3, reference implementation) — crossing = "cannot be
    arranged in a purely nested configuration"; nested pairs get order 0, knotted pairs order 1+.
  - **Wikipedia — Pseudoknot** (rank 4, cites Rivas & Eddy 1999; Antczak 2018) — qualitative definition
    ("half of one stem intercalated between the two halves of another"), "base pairing not well nested".

## Oracles and datasets

- **H-type `([)]` (positive):** pairs P1 = (0,2), P2 = (1,3). Crossing check `0 < 1 < 2 < 3` →
  **exactly one pseudoknot** whose two crossing pairs are P1, P2. (Antczak 2018 DBL example.)
- **Nested (negative control):** pairs (0,5) and (1,4): `0 < 1 < 4 < 5` fully nested → **no pseudoknot**.
- **Disjoint (negative control):** pairs (0,2) and (3,5): ranges `[0,2]` and `[3,5]` do not overlap →
  **no pseudoknot**.

## Invariants and edge cases

- **≥ 2 pairs required:** an empty pair set or a single pair can never cross → no pseudoknot (a derived
  consequence of the crossing definition, not an invented rule).
- **Endpoint normalization:** a pair whose endpoints are stored in either order (Position1 > Position2)
  is normalized via min/max to `open < close` **before** the crossing test — crossing is defined on
  `open < close` positions.
- **Determinism / order independence:** the same pair set always yields the same pseudoknots (a pure
  combinatorial O(n²) scan).
- **Property invariant (O(n²)):** for random pair sets, every reported pseudoknot's two pairs satisfy
  `i < k < j < l`.

## Deviations and assumptions

1. **Each crossing pair-of-pairs is reported as one Pseudoknot result** (crossing is a *binary relation*
   between two pairs). Grouping multiple mutually-crossing pairs into a single higher-**order**
   pseudoknot (Antczak 2018 DBL order assignment) is a richer feature that is **Not Implemented** and
   documented as such — a scope note, not an invented parameter.
2. **Empty / single-pair input returns no pseudoknots** — derived directly from "a pseudoknot needs two
   pairs that cross," not a free parameter.

**No source contradictions** — Antczak 2018, Smit & Knight 2008, biotite, and Wikipedia agree on the
crossing definition, the nested/disjoint negatives, and the symmetry of the crossing relation (which
pair is "the pseudoknot" is a labeling choice for removal/order; the *presence* of a crossing is
symmetric).
