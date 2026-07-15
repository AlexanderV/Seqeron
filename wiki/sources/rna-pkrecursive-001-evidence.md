---
type: source
title: "Evidence: RNA-PKRECURSIVE-001 (Recursive pknotsRG pseudoknot prediction — nested / multiple / over-arching knots)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md
sources:
  - docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md
source_commit: 562cd41454b50a8036d6440d5df044dff35a93c9
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-PKRECURSIVE-001

The validation-evidence artifact for test unit **RNA-PKRECURSIVE-001** — the **recursive-grammar
extension** of the canonical single-knot predictor: predict the energetically optimal RNA structure
that may contain **nested**, **multiple**, and **over-arching** pseudoknots, not just one. It is an
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the synthesizing concept is [[rna-pseudoknot-prediction]] (which this unit **enriches** — it is the
`PredictStructurePseudoknot` sibling `RnaSecondaryStructure` method that realizes the recursion the
single-knot unit explicitly did **not**). [[test-unit-registry]] tracks the unit.

This is the same **pknotsRG class** as [[rna-pkpredict-001-evidence|RNA-PKPREDICT-001]] — same energy
model, same 9.0 / 0.3 / 0.0 penalties, same canonization rules — but lifts the single-knot restriction:
it fills exactly the **PARTIAL-coverage gap** the single-knot evidence recorded.

## What "recursive" adds (the delta over RNA-PKPREDICT-001)

The pknotsRG grammar `knot = a ~~~ u ~~~ b ~~~ v ~~~ a' ~~~ w ~~~ b'` has three intervening loops
*u*, *v*, *w*. The single-knot predictor folds those loops **pseudoknot-free** (plain MFE). The
recursive class (Reeder & Giegerich 2004, verbatim) lets them "**fold internally in an arbitrary way,
including simple recursive pseudoknots**." The 2007 paper restates this and adds the whole-sequence
mechanism: the pseudoknot matrix value "**competes with values of unknotted foldings for the interval
(i, j)**" — per-interval competition is what lets the optimal structure carry **several** knots in
different regions and knots nested inside loops. Three realized capabilities:

1. **Recursion into loops** — a loop *u*/*v*/*w* may itself contain a further knot.
2. **Multiple knots** — the top level chains several knots along the sequence.
3. **Over-arching (enclosing) helix** — an outer helix may span a knot held in its loop.

## Sources (all as in RNA-PKPREDICT-001, re-used unchanged)

- **Reeder & Giegerich (2004)** *BMC Bioinformatics* 5:104 (rank 1, DOI 10.1186/1471-2105-5-104,
  PMC514697) — the recursive-class definition (loops fold internally including recursive knots),
  **O(n⁴) time / O(n²) space** (canonization reduces the independent boundaries from 8 → 4 versus
  Rivas & Eddy's O(n⁶)/O(n⁴)), the three canonization rules, and the penalties.
- **Reeder, Steffen & Giegerich (2007)** *Nucleic Acids Research* 35:W320 (rank 1, DOI
  10.1093/nar/gkm258, PMC1933184) — the per-interval **competition** with unknotted foldings and the
  restated recursive class (loops build secondary structure "including multiloops and pseudoknots").
- **pknotsRG reference source** `Energy.lhs` (rank 3, github.com/jensreeder/pknotsRG) — verbatim
  penalties **9.0 / 0.3 / 0.0**, same NN model for both helices, no extra per-pair penalty (unchanged
  from the single-knot unit).
- **Antczak et al. (2018)** *Bioinformatics* 34(8):1304 — the crossing condition `i < k < j < l`
  (`DetectPseudoknots`).

## Energy model and excluded classes

Energy is unchanged from the single-knot unit: **Turner NN stacking on both helices** + **9.0** knot
initiation + **0.3** per unpaired loop nt + **0.0** per pair inside the knot. **Excluded** from this
canonical simple-recursive class (verbatim): "more complex knotted structures like **triple crossing
helices or kissing hairpins** are excluded"; also bulged / unequal-length helices (canonization rule 1)
and chained / complex helix interactions.

## Oracles and datasets (all fully derivable)

- **Over-arching nested knot** — `AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` (38 nt): an A×8 5′ clamp +
  the designed H-type + a U×8 3′ clamp, so an 8-bp outer A·U helix over-arches the inner knot in its
  loop. Recursive structure `((((((((((((..[[[[..))))..]]]]))))))))`, **ΔG = −14.37 kcal/mol**,
  `HasPseudoknot == true`, 1 crossing. The single-knot `PredictStructurePseudoknot` and plain MFE both
  give −13.05 with `HasPseudoknot == false` (they cannot combine the outer helix **and** the inner
  knot). The over-arching knot is recovered **only** by the recursive method (−14.37 < −13.05).
- **Two separate knots** — an 80-nt sequence of two A·U-clamped copies of the designed H-type.
  Recursive structure recovers **two** crossing H-type knots (`DetectPseudoknots` crossing-count 32),
  **ΔG = −28.74**, `HasPseudoknot == true`; single-knot method and plain MFE both give −27.14 with
  **no** knot. Both knots recovered only by the recursive method.
- **No spurious knots** — plain hairpin `GGGGAAAACCCC` → recursive = MFE `((((....))))`, −5.28,
  `HasPseudoknot == false`, 0 crossings; A·U run `AUAUAUAUAUAUAUAU` → `((((((....))))))`, −0.26, no
  knot. Random sweep (seed 20260623, 150 seqs, len 12–38): recursive ΔG ≤ MFE ΔG for **all** (0
  violations), 0 spurious knots.
- **Single-knot parity** — `GGGGAACCCCAACCCCAAGGGG` (22 nt) → identical to
  `PredictStructurePseudoknot`: `((((..[[[[..))))..]]]]`, −8.76, same pairs. Recursion does not
  regress the single-knot case.

## Invariants

- **MFE fallback bound:** recursive `FreeEnergy ≤ CalculateMfeStructure(seq).FreeEnergy` for any
  sequence — the plain MFE is always in the search set (verified 0 violations on the random sweep).
- **No spurious knot:** the 9 kcal/mol initiation penalty gates a knot only when its two crossing
  helices beat both the penalty and the best pseudoknot-free alternative **for that interval**; plain
  inputs report `HasPseudoknot == false`.
- **Validity:** every index in range, each position paired ≤ once, reported knots have ≥1 genuine
  crossing (`i<k<j<l`).
- **Edge cases:** null / empty / too-short → empty pseudoknot-free structure (parity with the
  single-knot method); DNA spelling (T→U) folds identically; `minLoopSize < 3` clamps to 3.

## Deviations and assumptions

1. **PARTIAL coverage of the full pknotsRG recursion (documented, not an invented parameter).** The
   method realizes the recursive **class** — loops fold by the same recursive folder (a loop may
   contain a further knot), the top level chains multiple knots, and an enclosing helix may over-arch a
   knot. To stay tractable the H-type helices are enumerated by explicit start/end scan with maximal
   extension (rules 1–2) rather than the full 4-boundary ADP yield-parser; a knot is left-anchored
   within its interval before chaining; the enclosing-helix production is pursued only when the
   enclosed region is itself knotted (otherwise the pseudoknot-free Zuker MFE already covers it). The
   result is the faithful recursive class but **not** a guaranteed bit-identical reproduction of every
   yield the reference O(n⁴) parser would explore.
2. **Two-simultaneous-knot test cases are engineered, not random.** Two strong (G·C) H-type knots are
   the genuine MFE **only** when each knot region is isolated (flanking A·U clamps suppress the more
   stable cross-region nested alternative) — a property of the NN energy model, so the recovery is
   asserted on the engineered isolated-clamp sequence, not on random input.

**No source contradictions** — Reeder & Giegerich 2004/2007, the pknotsRG `Energy.lhs` source, and
Antczak 2018 agree on the recursive class, the 9.0 / 0.3 / 0.0 penalties, the same-NN-model helix
scoring, and the excluded complex classes. The only recorded items are the two PARTIAL / engineered
scope notes above.
