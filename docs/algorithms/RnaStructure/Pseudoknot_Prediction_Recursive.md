# RNA Pseudoknot Prediction — Recursive pknotsRG Grammar (nested / multiple knots)

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-PKRECURSIVE-001 (RNA-STRUCT-001 limitation fix) |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-23 |

## 1. Overview

Predicts an RNA secondary structure that may contain **multiple and recursively-nested canonical
H-type pseudoknots**, applying the recursive pknotsRG grammar of Reeder & Giegerich [1][2]
*throughout* the sequence rather than only to a single top-level knot. The whole sequence is folded
by a memoised interval recurrence in which a pseudoknot value "competes with values of unknotted
foldings for the interval (i, j)" [2], and the three loops of a knot may themselves contain further
pseudoknots [1]. It complements `PredictStructurePseudoknot` (single canonical H-type, left
unchanged) and reuses the identical energy model (Turner 2004 stacking + pknotsRG penalties 9.0 /
0.3 / 0.0 kcal/mol). It is a thermodynamic (energy-minimising) predictor for the canonical simple
*recursive* pseudoknot (csr-PK) class; it never returns a structure worse than the plain MFE.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A pseudoknot is a crossing arrangement of base pairs: pairs (i,j),(k,l) cross iff i<k<j<l [4].
Beyond a single H-type knot, biological RNAs may carry several pseudoknots, and a pseudoknot may sit
inside the loop of an outer helix. pknotsRG folds the entire sequence so the optimal structure can
contain such recursive arrangements within the canonical class [1][2].

### 2.2 Core Model

A *simple recursive pseudoknot* is two crossing helices with three loops, by the grammar [1]

```
knot = a ~~~ u ~~~ b ~~~ v ~~~ a' ~~~ w ~~~ b'
```

with the recursion rule, verbatim [1]: "we allow the unpaired strands u, v, w in a simple
pseudoknot to fold internally in an arbitrary way, **including simple recursive pseudoknots**."
The pknotsRG DP "extends the usual dynamic programming scheme for RNA folding" so that the
pseudoknot value "competes with values of unknotted foldings for the interval (i, j)" [2]; the
optimal structure over the whole sequence may therefore contain several knots and knots nested in
loops.

**Energy model (unchanged from the single-knot method).** ΔG(knot) = stacking(a·a') +
stacking(b·b') + ΔG(u) + ΔG(v) + ΔG(w) + Pᵢ + Pᵤₙₚ·(unpaired loop nucleotides), where the loop
energies ΔG(u/v/w) are obtained by the SAME recursive folder (so a loop may contain a further
knot), Pᵢ = 9.0 kcal/mol, Pᵤₙₚ = 0.3 kcal/mol, in-knot base pairs 0.0 [1][3].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-PKR-01 | Fixed temperature/buffer (Turner 2004, 37 °C) and pknotsRG penalties | Energies and knot/no-knot decisions shift with unmodelled conditions. |
| ASM-PKR-02 | Secondary-structure thermodynamics only; tertiary interactions not modelled | Tertiary-stabilised knots (e.g. BWYV) are not the predicted MFE. |
| ASM-PKR-03 | Two simultaneous strong knots are the MFE only when each region is isolated (e.g. A·U clamps) | On non-isolated G·C-rich inputs a cross-region nested fold can be more stable than two knots. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-PKR-01 | Returned ΔG ≤ ΔG of the plain MFE for the same sequence/parameters | The plain MFE is the fallback and is returned unchanged unless the recursive fold strictly improves it. |
| INV-PKR-02 | `HasPseudoknot` ⇒ the base-pair set contains ≥1 crossing pair (i<k<j<l) | Knot components are built from two helices that cross by construction. |
| INV-PKR-03 | Every position is paired at most once; all indices in [0, n) | Helices and recursively-folded sub-spans occupy disjoint, in-range index runs. |
| INV-PKR-04 | No spurious pseudoknot on a non-pseudoknotted sequence | A knot component is taken only when it lowers ΔG for its interval (9 kcal/mol penalty). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| rnaSequence | string | required | RNA or DNA (T read as U); case-insensitive | null/empty/too-short → empty result |
| minLoopSize | int | 3 | Minimum hairpin loop size for the nested folds | < 3 clamped to 3 (NNDB minimum) |

### 3.2 Output / Return Value

`PseudoknotStructure` record (shared with `PredictStructurePseudoknot`):

| Field | Type | Description |
|-------|------|-------------|
| Sequence | string | Folded sequence (upper-cased; T→U) |
| DotBracket | string | Two-layer dot-bracket: `()` = nested/outer helices, `[]` = crossing knot helices, `.` = unpaired |
| BasePairs | IReadOnlyList<(int,int)> | All base pairs as 0-based (5'<3') tuples, sorted by 5' position |
| FreeEnergy | double | ΔG° (kcal/mol) of the returned structure |
| HasPseudoknot | bool | True iff ≥1 crossing helix is present |

### 3.3 Preconditions and Validation

0-based indexing; upper-cased; T read as U. null/empty or shorter than the minimum canonical-knot
length (11 nt) returns an empty pseudoknot-free structure. Accepted alphabet A/C/G/U/T; other
characters do not pair.

## 4. Algorithm

### 4.1 High-Level Steps

1. Compute the plain pseudoknot-free MFE (`CalculateMfeStructure`) as baseline/fallback.
2. Fold the closed interval `[0, n−1]` with the memoised recurrence `F(i,j)`:
   - **Component 1 — pseudoknot-free block:** Zuker–Stiegler MFE of the sub-span (`CalculateMfeStructure`).
   - **Component 2 — H-type knot left-anchored at i:** scan the knot's 3' extent and inner
     boundaries (maximal helices, rules 1–2); score each loop u, v, w by `F` (recursive — a loop may
     contain a further knot); chain the remainder `F(r+1, j)`.
   - **Component 3 — over-arching helix:** pair (i,k) (maximal extension), fold the enclosed region
     `F(i+L, k−L)` recursively (so it can be knotted) and chain `F(k+1, j)`; pursued only when the
     enclosed region is itself knotted.
   - **Component 4 — leave i unpaired:** `F(i+1, j)`.
3. Take the minimum-energy decomposition per interval; memoise on (i,j).
4. Accept the recursive fold only if strictly below the plain MFE; otherwise return the plain MFE
   structure with `HasPseudoknot = false`. Render crossing helices on the `[]` layer.

### 4.2 Decision Rules / Reference Tables

| Constant | Value | Source |
|----------|-------|--------|
| Pseudoknot initiation Pᵢ | 9.0 kcal/mol | pknotsRG "creating a new pseudoknot: 9.0" [1][3] |
| Unpaired pseudoknot-loop base Pᵤₙₚ | 0.3 kcal/mol | pknotsRG "not paired base in pk: 0.3" [1][3] |
| Base pair inside pseudoknot | 0.0 kcal/mol | pknotsRG "basepair inside pseudoknot: 0.0" [3] |
| Helix stacking / terminal AU | Turner 2004 NN | reused `CalculateStemEnergy` (NNDB) |

No new energy parameter is introduced relative to `PredictStructurePseudoknot`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| PredictStructurePseudoknotRecursive | ~O(n⁴) (memoised intervals × per-interval helix scan; sub-span MFE folds) | O(n²) memo + O(n²) per MFE fold | Within the pknotsRG O(n⁴)/O(n²) envelope for the canonical class; intended for short-to-medium sequences. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.PredictStructurePseudoknotRecursive(string, int)`: public entry point.
- `RecursivePkFolder.Fold(i,j)`: memoised interval recurrence (Components 1–4).
- `RecursivePkFolder.TryKnotAnchoredAt(...)`, `EvaluateHTypeRecursive(...)`, `ScoreLoopRecursive(...)`:
  knot enumeration, knot scoring with recursive loops, and loop folding via `F`.

### 5.2 Current Behavior

Knot helices are scored by `CalculateStemEnergy` (Turner 2004). Each knot loop and each
enclosed/sub-span region folds by the same `Fold` recurrence, so pseudoknots may appear inside loops
(recursive class) and in multiple regions (top-level chain). Nucleotides left unpaired within a knot
loop are charged 0.3 kcal/mol. The recursive fold must beat the plain MFE by more than the DP
traceback tolerance to be accepted. No suffix tree is used: this is thermodynamic
scoring/enumeration, not exact-match search (see §5.3 of [Pseudoknot_Prediction](Pseudoknot_Prediction.md)).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Recursive class: knot loops u, v, w fold internally "in an arbitrary way, including simple
  recursive pseudoknots" — loops may contain further knots [1].
- Per-interval competition of pseudoknot vs unknotted foldings over the whole sequence, so the
  optimum can contain multiple knots and over-arching (nested) knots [2].
- Canonization rules 1 (equal-length, bulge-free helices) and 2 (maximal helix extent) [1].
- Pseudoknot initiation 9.0, unpaired-loop 0.3, in-knot base pair 0.0 kcal/mol [1][3].
- Turner 2004 nearest-neighbour stacking for all helices [3]; two-layer `()`/`[]` dot-bracket [4].

**Intentionally simplified:**

- Helix enumeration uses an explicit start/end scan with maximal extension and a left-anchored knot
  component before chaining, rather than the full 4-boundary ADP yield parser; **consequence:** the
  faithful recursive *class* (nested / multiple / over-arching knots) is produced, but not a
  guaranteed bit-identical reproduction of every yield the reference O(n⁴) parser explores.
- The over-arching-helix component is pursued only when the enclosed region is itself knotted;
  **consequence:** purely pseudoknot-free enclosing structure is left to the Zuker MFE (Component 1),
  which already covers it.
- Dangling/coaxial refinements at helix–loop junctions are not added separately; **consequence:** a
  predicted ΔG may differ slightly from pknotsRG's, but the knot/no-knot decision uses the same
  penalties and stacking model.

**Not implemented:**

- Kissing hairpins and triple-crossing / chained ("complex") helix interactions; **users should rely
  on:** pKiss / HotKnots / DotKnot — explicitly excluded from the canonical sr-PK class [1].
- Non-canonical (bulged, unequal-length) pseudoknot helices (canonization rule 1); **users should
  rely on:** the dedicated heuristics above.
- Tertiary-stabilised knots as the MFE; **users should rely on:** no current NN-thermodynamic
  alternative — this is an energy-model floor (ASM-PKR-02), not an algorithm gap.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Explicit helix scan vs full ADP yield parser | Assumption | Recursive class produced; not bit-identical to reference yields | accepted | §5.3; ASM-PKR-01. |
| 2 | Tertiary interactions not modelled | Assumption | Tertiary-stabilised knots not the MFE | accepted | ASM-PKR-02; inherent to NN predictors. |
| 3 | Two-knot optimum requires isolated regions | Assumption | Random G·C-rich inputs may favour a cross-region nested fold | accepted | ASM-PKR-03; engineered test cases. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | Empty structure (no pairs, ".", ΔG 0), HasPseudoknot=false | Contract parity with `CalculateMfeStructure`. |
| < 11 nt | No pseudoknot (returns plain MFE) | Shortest canonical knot = 11 nt. |
| Non-knot sequence / strong hairpin | HasPseudoknot=false; structure and ΔG = plain MFE | 9 kcal/mol penalty prevents spurious knots (INV-PKR-04). |
| Single canonical H-type | Identical result to `PredictStructurePseudoknot` | recursion does not regress the single-knot case. |
| DNA input (T) | Folded identically to the RNA spelling | T read as U. |

### 6.2 Limitations

Covers the canonical simple *recursive* pseudoknot (csr-PK) class: nested, multiple, and
over-arching canonical H-type knots. Out of scope: kissing hairpins, triple-crossing / chained
("complex") knots, and non-canonical (bulged / unequal-length) helices, all excluded from the
pknotsRG canonical class [1]. As a nearest-neighbour thermodynamic predictor it does not model
tertiary interactions, so tertiary-stabilised biological knots (e.g. BWYV, PDB 437D) are not
recovered as the MFE — an energy-model property of all NN-only predictors, not a defect.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Knot nested inside an outer (over-arching) helix:
var pk = RnaSecondaryStructure.PredictStructurePseudoknotRecursive("AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU");
// pk.HasPseudoknot == true
// pk.DotBracket    == "((((((((((((..[[[[..))))..]]]]))))))))"
// outer A·U helix (0,37)..(7,30); inner crossing knot stem1 (8,23)..(11,20), stem2 (14,29)..(17,26)
// pk.FreeEnergy (-14.37) < PredictStructurePseudoknot(...).FreeEnergy (-13.05)
```

### 7.2 Applications and Use Cases

- **Multi-pseudoknot RNAs:** viral 3' UTRs and some riboswitches carry more than one pseudoknot.
- **Over-arching knots:** a pseudoknot embedded within the loop of a larger hairpin/multiloop.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests.cs) — covers INV-PKR-01..04
- Evidence: [RNA-PKRECURSIVE-001-Evidence.md](../../../docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md)
- Related algorithms: [Pseudoknot_Prediction](Pseudoknot_Prediction.md) (single H-type), [Minimum_Free_Energy](Minimum_Free_Energy.md), [Pseudoknot_Detection](Pseudoknot_Detection.md)

## 8. References

1. Reeder J, Giegerich R. (2004). Design, implementation and evaluation of a practical pseudoknot folding algorithm based on thermodynamics. BMC Bioinformatics 5:104. https://doi.org/10.1186/1471-2105-5-104 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/)
2. Reeder J, Steffen P, Giegerich R. (2007). pknotsRG: RNA pseudoknot folding including near-optimal structures and sliding windows. Nucleic Acids Research 35:W320–W324. https://doi.org/10.1093/nar/gkm258 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1933184/)
3. Reeder J. pknotsRG source (Energy.lhs — pseudoknot penalties 9.0 / 0.3 / 0.0). https://github.com/jensreeder/pknotsRG
4. Antczak M, et al. (2018). New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. Bioinformatics 34(8):1304–1312. https://doi.org/10.1093/bioinformatics/btx783
