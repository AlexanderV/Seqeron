# RNA Secondary Structure Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-STRUCT-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

RNA secondary structure prediction identifies intramolecular base-pairing interactions and expresses the result as stems, loops, base-pair lists, and dot-bracket notation. In this repository, the public `PredictStructure` API is a heuristic structure builder rather than a full dynamic-programming fold: it enumerates candidate stem-loops, orders them by free energy, greedily selects a non-overlapping subset, and derives a dot-bracket string from the selected base pairs. The returned object also exposes a `Pseudoknots` field, but in the current public prediction path that field is effectively not meaningful because the greedy selection step excludes overlapping spans before pseudoknot detection runs. This makes the implementation useful for motif-oriented structure summaries while remaining simpler than exhaustive global RNA folding packages.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

RNA secondary structure is governed primarily by Watson-Crick and wobble base pairing. Common motifs include stems, hairpin loops, internal loops, bulges, multibranch loops, and pseudoknots. Dot-bracket notation is a standard textual representation in which `.` marks unpaired positions and paired positions are represented by matched brackets.

### 2.2 Core Model

The repository's structure predictor is built from three conceptual pieces confirmed in the source:

- Base-pair admissibility through canonical RNA pairing rules (`A-U`, `U-A`, `G-C`, `C-G`, and `G-U` / `U-G` wobble pairs).
- Candidate stem-loop enumeration over loop start and loop-size combinations.
- Greedy selection of the lowest-energy non-overlapping stem-loops to assemble the final structure.

The reported `MinimumFreeEnergy` in `SecondaryStructure` is the sum of the selected stem-loop energies rather than the direct output of the separate global `CalculateMinimumFreeEnergy` dynamic-programming routine.

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A useful structure summary can be formed from a non-overlapping subset of candidate stem-loops | The greedy selection may omit globally better combinations of interacting motifs |
| ASM-02 | Wobble pairing is biologically acceptable for the default prediction path | Disallowing wobble pairs would change the candidate set and predicted structures |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output `Sequence` is uppercase | `PredictStructure` uppercases the input sequence before further processing |
| INV-02 | `DotBracket.Length` equals `Sequence.Length` | The full dot-bracket generator allocates one character per sequence position |
| INV-03 | `MinimumFreeEnergy` equals the sum of `TotalFreeEnergy` over the selected stem-loops | The method computes the field by summing `selectedStemLoops.Sum(sl => sl.TotalFreeEnergy)` |
| INV-04 | Returned base pairs are ordered by ascending 5' position | The selected base pairs are materialized through `OrderBy(bp => bp.Position1)` |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Repository `PredictStructure` | Classical Nussinov / MFE folding |
|--------|-------------------------------|----------------------------------|
| Search strategy | Greedy selection over enumerated stem-loops | Dynamic programming over full sequence states |
| Optimization target | Non-overlapping low-energy stem-loops | Global optimality under the chosen recurrence and energy model |
| Pseudoknots | Detected after selection | Usually excluded from standard DP recurrences |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `rnaSequence` | `string` | required | RNA sequence to analyze | Empty input returns an empty `SecondaryStructure` |
| `minStemLength` | `int` | `3` | Minimum number of paired positions in a retained stem-loop | Applied during candidate filtering |
| `minLoopSize` | `int` | `3` | Minimum hairpin loop size | Passed to `FindStemLoops`, which clamps values below `3` |
| `maxLoopSize` | `int` | `10` | Maximum scanned hairpin loop size | Restricts candidate enumeration |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sequence` | `string` | Uppercased input sequence |
| `DotBracket` | `string` | Dot-bracket notation of the selected base pairs |
| `BasePairs` | `IReadOnlyList<BasePair>` | All selected base pairs ordered by position |
| `StemLoops` | `IReadOnlyList<StemLoop>` | Selected non-overlapping stem-loops |
| `Pseudoknots` | `IReadOnlyList<Pseudoknot>` | Crossing base-pair groups detected after selection |
| `MinimumFreeEnergy` | `double` | Sum of the selected stem-loop energies |

### 3.3 Preconditions and Validation

Empty input returns a `SecondaryStructure` with empty strings, empty collections, and zero energy. The method uppercases the sequence and does not perform separate alphabet validation. Stem-loop candidates are generated by `FindStemLoops` using default wobble pairing. If no candidate stem-loops survive selection, the method still returns a full-length dot-bracket string of `.` characters and zero total free energy.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return an empty structure when the input string is empty.
2. Uppercase the RNA sequence.
3. Enumerate candidate stem-loops with `FindStemLoops`.
4. Sort candidate stem-loops by `TotalFreeEnergy`.
5. Greedily select non-overlapping candidates in ascending-energy order.
6. Collect and order all base pairs from the selected stem-loops.
7. Generate full-length dot-bracket notation from the selected base pairs.
8. Sum the selected stem-loop energies and detect pseudoknots from the chosen base-pair list.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The predictor uses wobble pairing because `FindStemLoops` defaults `allowWobble` to `true` and `PredictStructure` does not override that default. Greedy selection is overlap-based: once a candidate stem-loop is selected, every sequence position between its `Start` and `End` is marked as used, and overlapping candidates are rejected. The generated dot-bracket output uses `(` and `)` for selected pairs and `.` for unpaired positions.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictStructure` | `O(C + h log h)` plus candidate search | `O(h + b)` | `C` = cost of `FindStemLoops`, `h` = number of candidate stem-loops, `b` = selected base pairs |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.PredictStructure(string, int, int, int)`: Builds the high-level structure summary.
- `RnaSecondaryStructure.FindStemLoops(string, int, int, int, bool)`: Supplies the candidate stem-loops used by the predictor.
- `RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int)`: Separate global MFE routine documented in [RNA_Free_Energy.md](./RNA_Free_Energy.md).

### 5.2 Current Behavior

The current `PredictStructure` implementation does not invoke the global dynamic-programming MFE routine when assembling the reported structure. Instead, it scores candidate stem-loops locally, sorts them by energy, and then performs a greedy non-overlap selection. Because that selection marks every position in a chosen stem-loop span as used, crossing stem combinations needed for pseudoknots are typically eliminated before the later `DetectPseudoknots(...)` call, so the returned `Pseudoknots` field is currently only a residual check over the surviving non-overlapping set. Although the general dot-bracket ecosystem includes additional bracket types for complex structures, the generated prediction string in this API uses only parentheses and dots.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- RNA Watson-Crick and wobble pairing rules through the base-pair helper logic.
- Dot-bracket output for the selected secondary-structure pairing pattern.

**Intentionally simplified:**

- Overall structure assembly uses greedy non-overlapping stem-loop selection rather than a global dynamic-programming optimum; **consequence:** the predicted structure is a heuristic summary, not a provably optimal fold.
- Pseudoknot detection runs only after greedy non-overlap selection; **consequence:** the public `PredictStructure` result currently does not preserve the crossing stem combinations needed for meaningful pseudoknot output, even though separate helper logic can detect crossings in arbitrary base-pair lists.

**Not implemented:**

- Comparative or consensus secondary-structure prediction; **users should rely on:** no current alternative in this class.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input sequence | Returns empty strings, empty collections, and zero energy | Explicit early return in `PredictStructure` |
| No candidate stem-loops | Dot-bracket contains only `.` characters and `MinimumFreeEnergy = 0` | No selected stem-loops contribute base pairs or energy |
| Overlapping low-energy candidates | Only the earliest accepted candidate in the energy-sorted greedy pass is retained | Overlap rejection is position-based |

### 6.2 Limitations

This predictor is a heuristic stem-loop assembler, not a full global RNA fold. It does not incorporate modified bases, comparative evidence, or pseudoknot-aware optimization, and its returned `Pseudoknots` field is effectively a byproduct of the already filtered non-overlapping structure rather than a full pseudoknot prediction. Because it depends on local stem-loop enumeration followed by greedy filtering, its output should be interpreted as a structured approximation rather than a full thermodynamic optimum.

## 7. Examples and Related Material

- [RNA-STRUCT-001](../../../tests/TestSpecs/RNA-STRUCT-001.md) documents the repository's RNA secondary-structure test specification.
- [RNA_Stemloop.md](./RNA_Stemloop.md) documents the stem-loop finder used as the candidate generator.
- [RNA_Free_Energy.md](./RNA_Free_Energy.md) documents the underlying energy calculations.

## 8. References

1. Nussinov, R., and A. B. Jacobson. 1980. Fast algorithm for predicting the secondary structure of single-stranded RNA. Proceedings of the National Academy of Sciences 77(11):6309-6313.
2. Zuker, M., and P. Stiegler. 1981. Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. Nucleic Acids Research 9(1):133-148.
3. Mathews, D. H., et al. 2004. Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proceedings of the National Academy of Sciences 101(19):7287-7292.
4. Wikipedia contributors. Nucleic acid secondary structure. Wikipedia.
5. Wikipedia contributors. Nussinov algorithm. Wikipedia.

