---
type: concept
title: "RNA pseudoknot prediction (canonical H-type, pknotsRG class)"
tags: [rna, algorithm]
sources:
  - docs/algorithms/RnaStructure/Pseudoknot_Prediction.md
  - docs/algorithms/RnaStructure/Pseudoknot_Prediction_Recursive.md
  - docs/Evidence/RNA-PKPREDICT-001-Evidence.md
  - docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md
  - docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md
source_commit: 2e3d94efc474888aead65cbf5c5ed6ea92d8ca6b
created: 2026-07-10
updated: 2026-07-16
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
**non-crossing** structures. Two test units realize this class: the **single**-knot predictor
`PredictStructurePseudoknot` (record [[rna-pkpredict-001-evidence]]) and its **recursive** extension
(record [[rna-pkrecursive-001-evidence]]) that composes **nested / multiple / over-arching** knots
in one structure (§6). [[test-unit-registry]] tracks the units and [[algorithm-validation-evidence]]
describes the artifact pattern.

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
  most once, and `DetectPseudoknots` finds ≥1 genuine crossing (`i<k<j<l`). That crossing check is the
  standalone [[rna-pseudoknot-detection|pseudoknot-detection]] primitive (RNA-PSEUDOKNOT-001) — a pure
  O(n²) scan over a given base-pair set, the analysis-side sibling of this energy-driven predictor.
- **Designed H-type oracle:** `GGGGAACCCCAACCCCAAGGGG` (22 nt) → `HasPseudoknot==true`, stem 1
  (0,15)(1,14)(2,13)(3,12) + stem 2 (6,21)(7,20)(8,19)(9,18), two-layer dot-bracket
  `((((..[[[[..))))..]]]]`, ΔG strictly below MFE.
- **Plain hairpin oracle:** `GGGGAAAACCCC` → `HasPseudoknot==false`, structure/ΔG equal MFE
  `((((....))))`.
- **Empty / null / too-short →** empty pseudoknot-free structure (no pairs, ΔG 0), parity with
  `CalculateMfeStructure`. DNA input (T→U) folds identically.

## 5. Scope, limitations, and relationships

A [[scientific-rigor|research-grade]] predictor. The **single**-knot unit
(`PredictStructurePseudoknot`, [[rna-pkpredict-001-evidence]]) predicts one canonical H-type knot; its
**recursive** sibling (§6, [[rna-pkrecursive-001-evidence]]) lifts that restriction. Two documented
boundaries remain (neither an invented parameter):

1. **PARTIAL pknotsRG coverage — split across the two units.** The full pknotsRG O(n⁴) grammar composes
   recursively-nested, over-arching, and multiple knots within one structure. The **single**-knot unit
   folds loops *u*/*v*/*w* with the pseudoknot-free MFE (`CalculateMinimumFreeEnergy`) and does **not**
   re-search a second knot inside a loop; the **recursive** unit (§6) does. Neither reproduces every
   yield of the reference 4-boundary ADP parser bit-for-bit (helices are enumerated by a maximal-extent
   start/end scan), but the recursive unit realizes the faithful recursive *class*.
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
Giegerich 2004/2007, the pknotsRG `Energy.lhs` source, the Wikipedia/Rivas & Eddy H-type geometry, the
BWYV structural record, and Antczak 2018's crossing condition are mutually consistent.

## 5a. Implementation surface (RNA-PKPREDICT-001 primary spec)

The **single**-knot unit is `RnaSecondaryStructure.PredictStructurePseudoknot(string rnaSequence,
int minLoopSize = 3)` in `Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs` (status **Simplified**);
`minLoopSize` is clamped up to the NNDB minimum of 3. It returns a **`PseudoknotStructure`** record
— `Sequence` (upper-cased, T→U), the two-layer `DotBracket`, `BasePairs` as 0-based `(5'<3')` tuples
sorted by 5' position, `FreeEnergy` (ΔG°, kcal/mol), and the `HasPseudoknot` flag. Internal helpers
`MaxHelixLength`, `EvaluateHType`, `ScoreLoop`, and `GeneratePseudoknotDotBracket` do maximal helix
extension, candidate scoring, loop folding, and two-layer rendering respectively.

The **enumeration** (spec §4.1) is not the reference ADP parser: it computes the plain
[[rna-minimum-free-energy-folding|MFE]] as baseline, then for each stem-1 start `i`/end maximal-extends
helix *a·a'* (canonization rules 1–2), and for each stem 1 enumerates stem-2 whose 5' strand starts
**inside loop 1** (between the two strands of *a*) and whose 3' strand ends **after** *a'* — so *b*
crosses *a* by construction (INV-PK-02). Each candidate is scored `stacking(a) + stacking(b) + 9.0 +
ΔG(u) + ΔG(v) + ΔG(w) + 0.3·(unpaired loop nt)`, where both stems reuse the module's
`CalculateStemEnergy` (Turner 2004 stacking + terminal AU/GU) and each loop span folds with
`CalculateMinimumFreeEnergy`/`CalculateMfeStructure`. The lowest-energy candidate is accepted only if
strictly below the plain MFE by more than the DP traceback tolerance. This makes the **implementation**
cost **O(n³) stem-start scan × loop-MFE, O(n²) space** — within, but tighter than, the pknotsRG
O(n⁴)/O(n²) envelope for the canonical class. No search data structure (a suffix tree was considered
and **not** used — this is thermodynamic enumeration, not exact-match search).

Two **intentional simplifications** distinguish the single-knot code from the reference: loops *u*/*v*/*w*
fold with the pseudoknot-free MFE rather than re-entering the full grammar (so a second knot nested
inside a loop is not predicted here — that is exactly what the **recursive** unit in §6 lifts), and
dangling/coaxial refinements at helix–loop junctions are not added separately (a predicted knot's ΔG may
differ slightly from pknotsRG's, but the knot/no-knot **decision** uses the same penalties and stacking
model). The shortest possible canonical knot is **11 nt** (2·2 stem pairs + 3 loop nt); shorter or
null/empty input returns the empty pseudoknot-free structure.

## 6. Recursive extension — nested / multiple / over-arching knots (RNA-PKRECURSIVE-001)

The **recursive** unit ([[rna-pkrecursive-001-evidence]]) realizes the part of the pknotsRG grammar the
single-knot unit leaves out. Same class, same energy model, same 9.0 / 0.3 / 0.0 penalties, same
canonization rules — but the three loops *u*/*v*/*w* now fold by the **same recursive folder** (a loop
may contain a further knot), the top level **chains multiple knots**, and an enclosing helix may
**over-arch** a knot in its loop. Reeder & Giegerich 2004 define the class verbatim: the unpaired
strands "fold internally in an arbitrary way, **including simple recursive pseudoknots**"; the 2007
paper supplies the whole-sequence mechanism — the pseudoknot value "**competes with values of unknotted
foldings for the interval (i, j)**," which is what lets one structure carry several knots and knots
nested inside loops.

Two constructed, fully-derivable oracles show the delta over the single-knot method (which recovers
neither, because it cannot combine an outer/second helix with a knot):

- **Over-arching knot** `AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` (38 nt, A×8 clamp + designed H-type +
  U×8 clamp) → `((((((((((((..[[[[..))))..]]]]))))))))`, **ΔG −14.37**, `HasPseudoknot==true`, 1
  crossing; single-knot method and plain MFE both give **−13.05** with no combined structure.
- **Two separate knots** (80 nt, two clamped H-type copies) → **two** crossing knots (crossing-count
  32), **ΔG −28.74**; single-knot method and MFE both **−27.14** with no knot.

The **MFE-fallback bound still holds** (recursive ΔG ≤ MFE; 0 violations on a 150-sequence random
sweep) and no spurious knots are reported on plain hairpins / A·U runs. Two engineered-scope notes
apply: two strong G·C knots are the genuine MFE **only** when each region is isolated by flanking A·U
clamps (else a cross-region nested helix is more stable — a property of the NN energy model), and the
implementation realizes the recursive *class* without guaranteeing every yield of the reference O(n⁴)
ADP parser. **Excluded** (verbatim): triple-crossing helices and kissing hairpins are not in the
canonical simple-recursive class.

### 6a. Implementation surface (RNA-PKRECURSIVE-001 primary spec)

The recursive unit is `RnaSecondaryStructure.PredictStructurePseudoknotRecursive(string rnaSequence,
int minLoopSize = 3)` in `Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs` (status **Simplified**;
`minLoopSize` clamped up to the NNDB minimum of 3). It returns the **same `PseudoknotStructure`
record** as the single-knot unit (§5a) — `Sequence` (upper-cased, T→U), the two-layer `DotBracket`
(`()` nested/over-arching helices, `[]` crossing knot helices), `BasePairs` as 0-based `(5'<3')`
tuples sorted by 5' position, `FreeEnergy` (ΔG°, kcal/mol), and `HasPseudoknot`. No new energy
parameter is introduced relative to `PredictStructurePseudoknot` (same 9.0 / 0.3 / 0.0 penalties, same
`CalculateStemEnergy` Turner 2004 stacking).

The **recursion contract** is a memoised interval recurrence `RecursivePkFolder.Fold(i,j)` over the
closed span `[0, n−1]`, memoised on `(i,j)` (so `HasPseudoknot`-per-interval competes with unknotted
foldings — the 2007-paper mechanism). Each interval takes the minimum-energy of **four competing
decompositions**:

1. **Pseudoknot-free block** — the Zuker–Stiegler [[rna-minimum-free-energy-folding|MFE]] of the
   sub-span (`CalculateMfeStructure`).
2. **H-type knot left-anchored at `i`** (`TryKnotAnchoredAt` → `EvaluateHTypeRecursive`) — scan the
   knot's 3' extent and inner boundaries by canonization rules 1–2 (maximal helices), score each loop
   *u*/*v*/*w* by `Fold` itself (**recursive** — a loop may hold a further knot, via
   `ScoreLoopRecursive`), then chain the remainder `Fold(r+1, j)`.
3. **Over-arching helix** — pair `(i,k)` at maximal extension, fold the enclosed region
   `Fold(i+L, k−L)` recursively (so it can be knotted) and chain `Fold(k+1, j)`; pursued **only when
   the enclosed region is itself knotted** (purely nested enclosures are already covered by
   Component 1, so no double-counting).
4. **Leave `i` unpaired** — `Fold(i+1, j)`.

The recursive fold is accepted only if strictly below the plain MFE by more than the DP-traceback
tolerance; otherwise the plain MFE is returned with `HasPseudoknot = false` (INV — no spurious knot).
The whole-sequence recurrence puts the **implementation** cost at **~O(n⁴)** (memoised intervals ×
per-interval helix scan × sub-span MFE folds), **O(n²)** memo + O(n²) per MFE fold — the full pknotsRG
O(n⁴)/O(n²) envelope, a step up from the single-knot unit's tighter O(n³) stem-scan (§5a) precisely
because loops and enclosures re-enter the folder. Intended for short-to-medium sequences. As in §5a,
no suffix tree is used (thermodynamic enumeration, not exact-match search). The single intentional
gap vs the reference: helices are enumerated by an explicit maximal-extent start/end scan with a
left-anchored knot component rather than the full 4-boundary ADP yield parser, so the recursive
*class* (nested / multiple / over-arching knots) is produced faithfully but not every reference yield
is reproduced bit-for-bit. Worked API oracle (spec §7.1): the over-arching-knot input
`AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` → `((((((((((((..[[[[..))))..]]]]))))))))`, ΔG **−14.37** <
`PredictStructurePseudoknot`'s **−13.05**.
