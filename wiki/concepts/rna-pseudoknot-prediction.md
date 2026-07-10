---
type: concept
title: "RNA pseudoknot prediction (canonical H-type, pknotsRG class)"
tags: [rna, algorithm]
sources:
  - docs/Evidence/RNA-PKPREDICT-001-Evidence.md
source_commit: b5cb44721acd570224bd0d138372f30e6f95b2a5
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-pkpredict-001-evidence
      evidence: "Test Unit ID: RNA-PKPREDICT-001 ... Algorithm: Pseudoknot Structure Prediction (canonical H-type, pknotsRG class)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-free-energy-turner-model
      source: rna-pkpredict-001-evidence
      evidence: "Energy model: Turner nearest-neighbour stacking for both pseudoknot helices, plus dangling/coaxial terms; pknotsRG Energy.lhs confirms helices use the SAME nearest-neighbour model as nested structures with no extra per-base-pair penalty (base pair inside pseudoknot: 0.0)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-minimum-free-energy-folding
      source: rna-pkpredict-001-evidence
      evidence: "Loops u/v/w fold with the existing pseudoknot-free MFE (CalculateMinimumFreeEnergy); invariant FreeEnergy ≤ CalculateMfeStructure(seq).FreeEnergy — the plain MFE is the always-available fallback baseline the predictor never underperforms"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-dot-bracket-notation
      source: rna-pkpredict-001-evidence
      evidence: "Designed canonical H-type sequence recovers the two-layer dot-bracket ((((..[[[[..))))..]]]] — crossing helices written with two independent bracket families"
      confidence: high
      status: current
---

# RNA pseudoknot prediction (canonical H-type, pknotsRG class)

The **crossing-helix layer** of the RNA secondary-structure family (test unit **RNA-PKPREDICT-001**):
predict the energetically optimal fold that may contain a single **pseudoknot** — two helices whose
base pairs **cross** (`i < k < j < l`, Antczak 2018). This is the one structural feature that every
**nested** predictor in the family is definitionally blind to: [[rna-minimum-free-energy-folding|MFE
folding]] and the [[rna-partition-function-mccaskill|McCaskill ensemble]] both enumerate only
**non-crossing** structures. The record is [[rna-pkpredict-001-evidence]]; [[test-unit-registry]]
tracks the unit and [[algorithm-validation-evidence]] describes the artifact pattern.

## 1. What a pseudoknot is (H-type geometry)

A **canonical simple recursive pseudoknot** (Reeder & Giegerich 2004) is "two crossing helices with
three intervening loops." In the standard **H-type** 5'→3' order (Rivas & Eddy 1999, via Wikipedia):

```
stem1-5' → loop1 → stem2-5' → loop2 → stem1-3' → loop3 → stem2-3'
```

Half of one stem is intercalated between the two halves of the other — that is what makes the pairs
cross. Written in [[rna-dot-bracket-notation|dot-bracket]] this needs a **two-layer** annotation,
using two independent bracket families: e.g. `((((..[[[[..))))..]]]]` — the `()` helix and the `[]`
helix interleave rather than nest. (A single-family notation cannot express a crossing; this is the
motivating example for WUSS's independent pairing systems.)

The pknotsRG grammar names the segments
`knot = a ~~~ u ~~~ b ~~~ v ~~~ a' ~~~ w ~~~ b'`: helix strand *a* pairs with *a'*, *b* with *b'*, and
*u*/*v*/*w* are the three separating loops (each may fold internally).

## 2. The canonical class (canonization rules)

pknotsRG restricts the search to a **canonical** subclass, which is exactly what makes the optimum
computable in **O(n⁴) time, O(n²) space** rather than the NP-hard general pseudoknot problem:

1. **Equal-length, bulge-free helices:** `|a| = |a'|`, `|b| = |b'|`; a helix is a contiguous run of
   pairs of equal strand length with no bulges.
2. **Maximal extent:** facing helices are extended to maximal Watson-Crick/GU length, so helix length
   is **determined by the sequence**, not searched independently.
3. **Fixed overlap boundary:** if two maximal helices would overlap, their boundary is fixed at an
   arbitrary point.

## 3. Energy model and penalties (sourced verbatim)

Both pseudoknot helices are scored with the **same** [[rna-free-energy-turner-model|Turner 2004
nearest-neighbour model]] as nested structures — the pknotsRG `Energy.lhs` reference source confirms
there is **no** extra per-base-pair penalty inside the knot. The pseudoknot-specific terms
(`Energy.lhs`, kcal/mol) are:

| Term | Value | Meaning |
|------|-------|---------|
| Pseudoknot initiation | **9.0** | "creating a new pseudoknot" — the destabilizing entry cost |
| Unpaired loop nucleotide | **0.3** | each unpaired nt inside a pseudoknot loop |
| Base pair inside pseudoknot | **0.0** | no extra penalty; scored by the NN model |

The **9.0 kcal/mol initiation penalty is the anti-spurious-knot gate:** a pseudoknot is reported only
when the two crossing helices' stabilization outweighs both this penalty **and** the best
pseudoknot-free alternative. Reeder & Giegerich chose 9.0 empirically ("performs better").

## 4. Invariants and oracles

- **INV (MFE fallback bound):** `FreeEnergy ≤ CalculateMfeStructure(seq).FreeEnergy` — the plain
  [[rna-minimum-free-energy-folding|MFE fold]] is always in the search set, so the predictor never
  returns a worse structure.
- **INV (no spurious knot):** when a nested structure is more stable, `HasPseudoknot == false` and the
  output equals the plain MFE.
- **INV (validity):** when a knot is returned, every index is in range, each position is paired at
  most once, and `DetectPseudoknots` finds ≥1 genuine crossing (`i<k<j<l`).
- **Designed H-type oracle:** `GGGGAACCCCAACCCCAAGGGG` (22 nt) → `HasPseudoknot==true`, stem 1
  (0,15)(1,14)(2,13)(3,12) + stem 2 (6,21)(7,20)(8,19)(9,18), two-layer dot-bracket
  `((((..[[[[..))))..]]]]`, ΔG strictly below MFE.
- **Plain hairpin oracle:** `GGGGAAAACCCC` → `HasPseudoknot==false`, structure/ΔG equal MFE
  `((((....))))`.
- **Empty / null / too-short →** empty pseudoknot-free structure (no pairs, ΔG 0), parity with
  `CalculateMfeStructure`. DNA input (T→U) folds identically.

## 5. Scope, limitations, and relationships

A [[scientific-rigor|research-grade]] predictor of the **single canonical H-type** pseudoknot. Two
documented boundaries (neither an invented parameter):

1. **PARTIAL pknotsRG coverage.** The full pknotsRG O(n⁴) grammar also composes recursively-nested and
   over-arching / multiple knots within one structure; only the **single** canonical H-type is
   implemented. Loops *u*/*v*/*w* fold with the existing pseudoknot-free MFE
   (`CalculateMinimumFreeEnergy`) rather than re-searching a second knot inside a loop.
2. **Tertiary-stabilized knots are not recovered.** Real knots such as the **BWYV** frameshifting
   pseudoknot (PDB 437D, `GGCGCGGCACCGUCCGCGGAACAAACGG`) are held together by minor-groove triplexes
   and Mg²⁺/Na⁺ coordination — interactions **outside** the secondary-structure NN model. The NN
   optimum here is the pseudoknot-free stem-1 hairpin, not the crystallographic knot: a documented
   limit of **all** NN-only pseudoknot predictors, recorded to prevent over-claiming.

Within the family this unit **depends on** [[rna-free-energy-turner-model]] (both helices are scored
with its Turner terms) and on [[rna-minimum-free-energy-folding]] (the plain MFE is the fallback
baseline and folds the internal loops); it emits a two-layer [[rna-dot-bracket-notation]] structure
and builds on the [[rna-base-pairing]] Watson-Crick + G-U wobble rule. It is the crossing-structure
complement of the nested [[rna-minimum-free-energy-folding|MFE folder]] and
[[rna-partition-function-mccaskill|McCaskill ensemble]]. **No source contradictions** — Reeder &
Giegerich 2004, the pknotsRG `Energy.lhs` source, the Wikipedia/Rivas & Eddy H-type geometry, the
BWYV structural record, and Antczak 2018's crossing condition are mutually consistent.
