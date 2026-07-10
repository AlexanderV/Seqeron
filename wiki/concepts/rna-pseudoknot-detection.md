---
type: concept
title: "RNA pseudoknot detection (crossing base pairs in a given structure)"
tags: [rna, algorithm]
sources:
  - docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md
source_commit: ae0dfc54b6b719a8fa68c2f120f3f4e3235cd02e
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-pseudoknot-001-evidence
      evidence: "Test Unit ID: RNA-PSEUDOKNOT-001 ... Algorithm: Pseudoknot Detection (crossing base pairs)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-pseudoknot-prediction
      source: rna-pseudoknot-001-evidence
      evidence: "Detection is the crossing primitive the predictor's validity invariant leans on (DetectPseudoknots finds ≥1 genuine crossing i<k<j<l when a knot is returned); both facets share the Antczak 2018 crossing condition, but detection scans a given base-pair set while prediction folds a sequence."
      confidence: high
      status: current
---

# RNA pseudoknot detection (crossing base pairs in a given structure)

The **detection / analysis** facet of the RNA crossing-helix family (test unit
**RNA-PSEUDOKNOT-001**, `RnaSecondaryStructure.DetectPseudoknots`): given a **set of base pairs** that
already describes a structure, find the **crossing** pairs that make it a **pseudoknot**. This is a
pure **O(n²) combinatorial scan** — no sequence, no energy model — and is therefore genuinely distinct
from its sibling [[rna-pseudoknot-prediction|pknotsRG predictor]], which takes a *sequence* and *folds*
it (O(n⁴), Turner nearest-neighbour energy). Detection is the primitive the predictor's validity
invariant relies on. [[test-unit-registry]] tracks the unit; [[algorithm-validation-evidence]] describes
the artifact pattern; the source record is [[rna-pseudoknot-001-evidence]].

## 1. The crossing condition (the definitional core)

Two base pairs (i,j) and (k,l), each written **open < close**, **cross** exactly when

```
i < k < j < l
```

— one pair's opener lies inside the other pair while its closer lies outside (Antczak et al. 2018,
verbatim: "for any pair (*i, i'*) there exists another one (*j, j'*) such that *i < j < i' < j'*").
Crossing is what "not well nested" means: base pairs "overlap one another in sequence position"
(Wikipedia / Rivas & Eddy 1999). The two exhaustive negatives are:

| Relation | Condition | Pseudoknot? |
|----------|-----------|-------------|
| **Crossing** | `i < k < j < l` | **yes** |
| **Nested** | `i < k < l < j` (one pair fully inside the other) | no |
| **Disjoint** | `j < k` (side-by-side, non-overlapping ranges) | no |

Written in [[rna-dot-bracket-notation|dot-bracket]] the minimal H-type crossing is `([)]` — the `()`
and `[]` families interleave rather than nest, which is precisely why the notation needs independent
bracket families per crossing layer.

## 2. Contract, invariants, and edge cases

- **Input is a base-pair set, output is the crossings.** Each **crossing pair-of-pairs is reported as
  one pseudoknot** — crossing is a *binary relation* between two pairs.
- **≥ 2 pairs required:** an empty set or a single pair can never cross → no pseudoknot (a derived
  consequence of the crossing definition, not an invented rule).
- **Endpoint normalization:** a pair stored in either order (`Position1 > Position2`) is normalized via
  min/max to `open < close` **before** the crossing test.
- **Determinism / order independence:** the same pair set always yields the same pseudoknots.
- **Property invariant (O(n²)):** for any random pair set, every reported pseudoknot's two pairs satisfy
  `i < k < j < l`.

Worked oracles (Antczak 2018 / Wikipedia): `([)]` = (0,2)+(1,3) → `0<1<2<3` → **one** pseudoknot;
nested (0,5)+(1,4) → `0<1<4<5` → **none**; disjoint (0,2)+(3,5) → ranges do not overlap → **none**.

## 3. Scope, limitations, and relationships

A [[scientific-rigor|research-grade]] combinatorial detector. One documented boundary (a scope note,
not an invented parameter): **pseudoknot-order assignment is Not Implemented.** Grouping multiple
mutually-crossing pairs into a single higher-**order** pseudoknot (Antczak 2018's DBL order layering —
order 0 `()`, 1 `[]`, 2 `{}`, 3 `<>`, 4–8 letters; Smit et al. 2008 pseudoknot-removal / order
assignment) is a richer feature this unit does not compute; it reports the crossing *relations*
themselves.

Within the family, detection is the crossing-check primitive behind the
[[rna-pseudoknot-prediction|predictor]]'s validity invariant (when the predictor returns a knot,
`DetectPseudoknots` confirms ≥ 1 genuine crossing) — the two facets **share** the Antczak 2018 crossing
condition `i < k < j < l` but differ completely in input (base-pair set vs sequence), machinery
(combinatorial scan vs energy DP) and cost (O(n²) vs O(n⁴)). It operates on the pairs produced by the
[[rna-dot-bracket-notation]] parser and complements the nested-only
[[rna-minimum-free-energy-folding|MFE folder]] / [[rna-partition-function-mccaskill|McCaskill ensemble]],
which are definitionally blind to crossings. **No source contradictions** — Antczak 2018, Smit & Knight
2008, biotite, and Wikipedia agree on the crossing definition, the nested/disjoint negatives, and the
symmetry of the crossing relation.
