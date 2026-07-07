# RNA Pseudoknot Structure Prediction (canonical H-type, pknotsRG class)

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-PKPREDICT-001 (RNA-STRUCT-001 limitation fix) |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-23 |

## 1. Overview

Predicts an RNA secondary structure that may contain a single **canonical H-type pseudoknot** —
two crossing helices with three intervening loops — using the *canonical simple recursive
pseudoknot* class of Reeder & Giegerich's pknotsRG algorithm [1]. The two crossing helices are
scored with the same Turner 2004 nearest-neighbour stacking model used by the pseudoknot-free MFE
folder in this module; the three connecting loops fold independently with that MFE; and the
pknotsRG pseudoknot-specific penalties are added (initiation 9.0 kcal/mol; 0.3 kcal/mol per
unpaired loop nucleotide) [1][2]. The candidate H-type fold is accepted only when its total free
energy is strictly below the plain pseudoknot-free MFE, so no spurious pseudoknots are introduced.
It is a thermodynamic (energy-minimising) predictor, exact for the canonical single-knot class it
covers and heuristic-free within that class.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A pseudoknot is a non-nested ("crossing") arrangement of base pairs: pairs (i,j) and (k,l) cross
iff i < k < j < l [5]. The most common biological case is the **H-type** pseudoknot, where the
unpaired loop of a hairpin base-pairs with a single-stranded region outside the hairpin, so half of
one helix is intercalated between the two halves of the other [3]. Standard O(n³) nearest-neighbour
folders (Zuker–Stiegler) are pseudoknot-free by construction; predicting pseudoknots requires a
different algorithm class [1].

### 2.2 Core Model

pknotsRG defines a *simple recursive pseudoknot* as two crossing helices with three intervening
loops, written by the grammar [1]

```
knot = a ~~~ u ~~~ b ~~~ v ~~~ a' ~~~ w ~~~ b'
```

read 5'→3', where helix strand `a` pairs with `a'`, helix strand `b` pairs with `b'`, and `u`,
`v`, `w` are the three loops. The H-type ordering is therefore
**stem1-5' (a) · loop1 (u) · stem2-5' (b) · loop2 (v) · stem1-3' (a') · loop3 (w) · stem2-3' (b')**
[3]. Because `b` lies between the two strands of `a` while `b'` lies outside `a'`, the two helices
cross — a genuine pseudoknot.

**Energy model.** ΔG(knot) = stacking(a·a') + stacking(b·b') + ΔG(u) + ΔG(v) + ΔG(w) + Pᵢ +
Pᵤₙₚ·(unpaired loop nucleotides), where stacking(·) is the Turner 2004 nearest-neighbour stack/
terminal-AU energy [reused from the module's `CalculateStemEnergy`], ΔG(loop) is the pseudoknot-free
MFE of the loop span, Pᵢ = 9.0 kcal/mol is the pseudoknot initiation penalty and Pᵤₙₚ = 0.3
kcal/mol is the penalty per unpaired loop nucleotide [1][2]. Base pairs inside the pseudoknot carry
no extra term ("basepair inside pseudoknot: 0.0") [2].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-PK-01 | Fixed temperature/buffer (Turner 2004, 37 °C); pseudoknot penalties as published by pknotsRG | Energies and the knot/no-knot decision shift with conditions not modelled (e.g. ionic strength). |
| ASM-PK-02 | Secondary-structure thermodynamics only; tertiary interactions (triplexes, Mg²⁺ coordination) are not modelled | Tertiary-stabilised knots (e.g. BWYV) are not the predicted MFE structure. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-PK-01 | Returned ΔG ≤ ΔG of the plain MFE structure for the same sequence/parameters | The plain MFE is the baseline and is returned unchanged unless a candidate knot strictly improves it. |
| INV-PK-02 | `HasPseudoknot` ⇒ the base-pair set contains a crossing pair (i<k<j<l) | The accepted candidate is built from two helices that cross by construction (stem2-5' inside stem1, stem2-3' outside). |
| INV-PK-03 | Every position is paired at most once; all indices in [0, n) | Helices use disjoint, in-range index runs and loop folds operate on disjoint sub-spans. |
| INV-PK-04 | No spurious pseudoknot | A knot is accepted only when its energy is strictly below the MFE (by more than the DP tolerance). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| rnaSequence | string | required | RNA or DNA sequence (T read as U); case-insensitive | null/empty/too-short → empty result |
| minLoopSize | int | 3 | Minimum hairpin loop size for loop folding | < 3 is clamped to 3 (NNDB minimum) |

### 3.2 Output / Return Value

`PseudoknotStructure` record:

| Field | Type | Description |
|-------|------|-------------|
| Sequence | string | Folded sequence (upper-cased; T→U) |
| DotBracket | string | Two-layer dot-bracket: `()` = stem 1, `[]` = crossing stem 2, `.` = unpaired; single-family MFE notation when no knot |
| BasePairs | IReadOnlyList<(int,int)> | All base pairs as 0-based (5'<3') tuples, sorted by 5' position |
| FreeEnergy | double | ΔG° (kcal/mol) of the returned structure |
| HasPseudoknot | bool | True iff the returned structure contains a crossing helix |

### 3.3 Preconditions and Validation

0-based indexing. Sequence is upper-cased and T is read as U (A–U/A–T identical for folding).
null/empty or shorter than the minimum canonical-knot length (11 nt) returns an empty
pseudoknot-free structure (no pairs, all dots, ΔG = 0). Accepted alphabet A/C/G/U/T; other
characters simply do not pair.

## 4. Algorithm

### 4.1 High-Level Steps

1. Compute the plain pseudoknot-free MFE structure (`CalculateMfeStructure`) as baseline/fallback.
2. Enumerate stem-1 (a·a'): for each 5' start i and each 3' end, maximal-extend the helix
   (canonization rules 1–2) to length L1.
3. For each stem 1, enumerate stem-2 (b·b'): its 5' strand starts in loop 1 (between the two
   strands of a) and its 3' strand ends after a', maximal-extended to length L2, so b crosses a.
4. Score the candidate: stacking(a) + stacking(b) + 9.0 + fold each of the three loops with the
   pseudoknot-free MFE + 0.3 per unpaired loop nucleotide.
5. Keep the lowest-energy candidate; accept it only if strictly below the plain MFE; otherwise
   return the plain MFE structure with `HasPseudoknot = false`.
6. Render the accepted knot in two-layer dot-bracket (`()` stem 1, `[]` stem 2).

### 4.2 Decision Rules / Reference Tables

| Constant | Value | Source |
|----------|-------|--------|
| Pseudoknot initiation Pᵢ | 9.0 kcal/mol | pknotsRG: "creating a new pseudoknot: 9.0" [1][2] |
| Unpaired pseudoknot-loop base Pᵤₙₚ | 0.3 kcal/mol | pknotsRG: "not paired base in pk: 0.3" [1][2] |
| Base pair inside pseudoknot | 0.0 kcal/mol | pknotsRG: "basepair inside pseudoknot: 0.0" [2] |
| Helix stacking / terminal AU | Turner 2004 NN | reused `CalculateStemEnergy` (NNDB) |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| PredictStructurePseudoknot | O(n³) stem-start scan × loop MFE | O(n²) | Within the pknotsRG O(n⁴)/O(n²) envelope for the canonical single-knot class; loop folds reuse the module MFE. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.PredictStructurePseudoknot(string, int)`: returns `PseudoknotStructure`; the public entry point for this unit.
- `MaxHelixLength(...)`, `EvaluateHType(...)`, `ScoreLoop(...)`, `GeneratePseudoknotDotBracket(...)`: internal helpers (maximal helix extension, candidate scoring, loop folding, two-layer rendering).

### 5.2 Current Behavior

The two stems are scored by passing their base-pair lists to the existing `CalculateStemEnergy`
(Turner 2004 stacking + terminal AU/GU). Each loop span folds with `CalculateMinimumFreeEnergy`/
`CalculateMfeStructure`; nucleotides left unpaired by that fold are charged 0.3 kcal/mol. The
candidate must beat the plain MFE by more than the DP traceback tolerance to be accepted. No search
data structure (suffix tree) is used: this is a thermodynamic scoring/enumeration, not exact-match
search — see §5.3.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Canonical H-type geometry: two crossing helices a·a', b·b' with three loops u, v, w [1][3].
- Canonization rule 1 (equal-length, bulge-free helices) and rule 2 (maximal helix extent) [1].
- Pseudoknot initiation 9.0, unpaired-loop 0.3, in-knot base pair 0.0 kcal/mol [1][2].
- Turner 2004 nearest-neighbour stacking for both helices [2].
- Two-layer `()`/`[]` dot-bracket for crossing helices (ViennaRNA/WUSS) [5].

**Intentionally simplified:**

- Loop folding: loops fold with the pseudoknot-free MFE rather than re-entering the full pknotsRG
  grammar; **consequence:** a second pseudoknot nested inside a loop is not predicted.
- Dangling/coaxial refinements at the helix–loop junctions are not added separately;
  **consequence:** ΔG of a predicted knot may differ slightly from pknotsRG's, but the knot/no-knot
  decision uses the same penalties and stacking model.

**Not implemented:**

- The full pknotsRG O(n⁴) recursion over recursively-nested / multiple / over-arching pseudoknots
  in one structure; **users should rely on:** pknotsRG / pknots / IPknot for those classes. This
  unit covers the single canonical H-type knot only.
- Kissing-hairpin and non-canonical (bulged, unequal-length) pseudoknots; **users should rely on:**
  the dedicated heuristics (HotKnots, DotKnot) — out of scope of the canonical class.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Single canonical H-type knot only | Assumption | No recursive/multiple knots predicted | accepted | Documented PARTIAL; see ASM-PK in Evidence and §6.2. |
| 2 | Tertiary interactions not modelled | Assumption | Tertiary-stabilised knots (BWYV) not the MFE | accepted | Inherent to all NN-only predictors (ASM-PK-02). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | Empty structure (no pairs, ".", ΔG 0), HasPseudoknot=false | Contract parity with `CalculateMfeStructure`. |
| < 11 nt | No pseudoknot (returns plain MFE) | Shortest canonical knot = 2·2 + 2·2 + 3 loops = 11 nt. |
| Non-knot sequence / strong hairpin | HasPseudoknot=false; structure and ΔG = plain MFE | 9 kcal/mol penalty prevents spurious knots (INV-PK-04). |
| DNA input (T) | Folded identically to the RNA spelling | T read as U. |

### 6.2 Limitations

Covers the **single canonical H-type pseudoknot** (pknotsRG canonical simple recursive class).
Recursively-nested pseudoknots, multiple/over-arching knots in one structure, kissing hairpins,
and non-canonical (bulged / unequal-length) helices are out of scope. As a nearest-neighbour
thermodynamic predictor it does not model tertiary interactions, so tertiary-stabilised biological
knots (e.g. the BWYV frameshifting pseudoknot, PDB 437D) are not recovered as the MFE structure —
a documented property of NN-only pseudoknot prediction, not a defect.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var pk = RnaSecondaryStructure.PredictStructurePseudoknot("GGGGAACCCCAACCCCAAGGGG");
// pk.HasPseudoknot == true
// pk.DotBracket    == "((((..[[[[..))))..]]]]"
// stem 1: (0,15)(1,14)(2,13)(3,12); stem 2 (crossing): (6,21)(7,20)(8,19)(9,18)
// pk.FreeEnergy    <  CalculateMfeStructure("GGGGAACCCCAACCCCAAGGGG").FreeEnergy
```

### 7.2 Applications and Use Cases

- **Frameshifting / readthrough signals:** H-type pseudoknots stimulate −1 ribosomal frameshifting in many viruses.
- **Telomerase / catalytic RNAs:** the telomerase P2b-P3 pseudoknot is an H-type knot.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs) — covers INV-PK-01..04
- Evidence: [RNA-PKPREDICT-001-Evidence.md](../../../docs/Evidence/RNA-PKPREDICT-001-Evidence.md)
- Related algorithms: [Minimum_Free_Energy](Minimum_Free_Energy.md), [Pseudoknot_Detection](Pseudoknot_Detection.md)

## 8. References

1. Reeder J, Giegerich R. (2004). Design, implementation and evaluation of a practical pseudoknot folding algorithm based on thermodynamics. BMC Bioinformatics 5:104. https://doi.org/10.1186/1471-2105-5-104
2. Reeder J. pknotsRG source (Energy.lhs — pseudoknot penalties 9.0 / 0.3 / 0.0). https://github.com/jensreeder/pknotsRG
3. Pseudoknot. Wikipedia (cites Rivas & Eddy 1999; Staple & Butcher 2005). https://en.wikipedia.org/wiki/Pseudoknot
4. Su L, Chen L, Egli M, Berger JM, Rich A. (1999). Minor groove RNA triplex in the crystal structure of a ribosomal frameshifting viral pseudoknot. Nat Struct Biol 6(3):285–292. https://pmc.ncbi.nlm.nih.gov/articles/PMC7097825/ (PDB 437D: https://www.rcsb.org/structure/437D)
5. Antczak M, et al. (2018). New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. Bioinformatics 34(8):1304–1312. https://doi.org/10.1093/bioinformatics/btx783
