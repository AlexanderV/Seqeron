# RNA Secondary Structure Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-STRUCT-001 |
| Related Projects | N/A |
| Implementation Status | Production |
| Last Reviewed | 2026-06-23 |

## 1. Overview

RNA secondary structure prediction identifies intramolecular base-pairing interactions and expresses the result as stems, loops, base-pair lists, and dot-bracket notation. This repository exposes **two** prediction paths over the same Turner 2004 energy model:

- **MFE-optimal (default for correctness):** `CalculateMfeStructure` / `PredictStructureMfe` return the globally optimal pseudoknot-free structure by Zuker–Stiegler (1981) dynamic-programming **traceback** over the same V/W/WM matrices used to compute the scalar MFE. The reconstructed structure's free energy equals the value returned by `CalculateMinimumFreeEnergy` for the same input — they are mutually consistent by construction [2].
- **Stem-loop assembly (restricted):** `PredictStructure` enumerates hairpin stem-loop candidates (`FindStemLoops`: perfect stems ≥ `minStemLength`, loops `minLoopSize`..`maxLoopSize`) and selects the minimum-free-energy non-overlapping subset exactly, by the Zuker–Stiegler exterior-loop recursion over the candidates (RNA-MFE-001, 2026-11-25; previously a greedy pass that could select destabilising stem-loops). Its reported energy is the free energy of the returned structure (`CalculateFreeEnergy`, ≤ 0), but the search space has no multiloops, bulges or interior loops, so it is generally **not** the global optimum.

Previously only the stem-loop path was available, so the returned structure was an approximation even though the scalar MFE was correct. The traceback closes that gap (MCP `rna_mfe_structure`).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

RNA secondary structure is governed primarily by Watson-Crick and wobble base pairing. Common motifs include stems, hairpin loops, internal loops, bulges, multibranch loops, and pseudoknots. Dot-bracket notation is a standard textual representation in which `.` marks unpaired positions and paired positions are represented by matched brackets.

### 2.2 Core Model

**MFE-optimal path (`CalculateMfeStructure`).** The Zuker–Stiegler recurrences fill three matrices [2] (taught with explicit equations in MIT 6.047 Lecture 08, Fig. 13 — the F/C/M/M¹ decomposition [6]):

- `V(i,j)` — minimum energy of `[i..j]` given `i·j` paired: the minimum over a hairpin, a stack `V(i+1,j-1)`, an interior/bulge loop `V(i',j')`, and a multiloop `WM(i+1,j-1)`.
- `WM(i,j)` — energy of a multiloop region with ≥ 1 branch (a branch `V(i,j)` + its d2 stem term, an unpaired end, or a split); a multiloop closed by `(i,j)` needs ≥ 2 branches: `Mclose(i,j) + min_k WM(i+1,k) + WM(k+1,j−1)`.
- Interior loops use the Turner 2004 int11/int21/int22 tables and 1×n / 2×3 / generic mismatch rules exactly as ViennaRNA 2.7.2 `vrna_E_internal`; exterior and multiloop stems use the d2 dangle model (terminal mismatch when both neighbours exist). With NNDB-precision parameters the MFE and structure equal ViennaRNA `RNA.fold` (RNA-STRUCT-001 Evidence, Validation 2026-11).
- `W(j)` — overall MFE of the prefix `[0..j]`: `j` unpaired, or a helix `V(k,j)` with `W(k-1)` before it.

The **traceback** re-derives, at each cell, which recurrence option attained the stored optimum and recurses into the corresponding sub-problem(s), recording a base pair whenever a `V(i,j)` cell is entered. Because it mirrors the fill exactly, the reconstructed structure's energy equals `W(0,n-1)` — the scalar MFE.

**Stem-loop assembly path (`PredictStructure`).** Built from three conceptual pieces:

- Base-pair admissibility through canonical RNA pairing rules (`A-U`, `U-A`, `G-C`, `C-G`, and `G-U` / `U-G` wobble pairs).
- Candidate stem-loop enumeration over loop start and loop-size combinations.
- Exact selection of the lowest-free-energy set of non-overlapping stem-loops: `W(j) = min{W(j−1), min_c W(c.Start−1) + E(c)}`, E(c) = free energy of the sequence folded into candidate c alone (under d2 the exterior terms depend only on the sequence, so disjoint candidates add).

The reported `MinimumFreeEnergy` of `PredictStructure` is `CalculateFreeEnergy` of the returned structure (≥ the global `CalculateMinimumFreeEnergy`, ≤ 0).

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A useful structure summary can be formed from a non-overlapping subset of candidate stem-loops | Structures needing multiloops / interior loops / nested hairpins are outside the search space |
| ASM-02 | Wobble pairing is biologically acceptable for the default prediction path | Disallowing wobble pairs would change the candidate set and predicted structures |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output `Sequence` is uppercase | `PredictStructure` uppercases the input sequence before further processing |
| INV-02 | `DotBracket.Length` equals `Sequence.Length` | The full dot-bracket generator allocates one character per sequence position |
| INV-03 | `MinimumFreeEnergy` equals the sum of `TotalFreeEnergy` over the selected stem-loops | The method computes the field by summing `selectedStemLoops.Sum(sl => sl.TotalFreeEnergy)` |
| INV-04 | Returned base pairs are ordered by ascending 5' position | The selected base pairs are materialized through `OrderBy(bp => bp.Position1)` (stem-loop assembly) / sorted by `Position1` (traceback) |
| INV-05 | `CalculateMfeStructure(s).FreeEnergy == CalculateMinimumFreeEnergy(s)` | The traceback follows the same V/W/WM recurrences that produced the scalar MFE |
| INV-06 | The traceback structure is pseudoknot-free (no crossing pairs) | The Zuker–Stiegler recurrences only nest or juxtapose helices |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Repository `PredictStructure` | Classical Nussinov / MFE folding |
|--------|-------------------------------|----------------------------------|
| Search strategy | Exact exterior-loop recursion over enumerated stem-loops | Dynamic programming over full sequence states |
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

The predictor uses wobble pairing because `FindStemLoops` defaults `allowWobble` to `true` and `PredictStructure` does not override that default. Selection is over intervals: chosen stem-loops occupy disjoint `[Start, End]` spans (no nesting), so overlapping candidates are never combined. The generated dot-bracket output uses `(` and `)` for selected pairs and `.` for unpaired positions.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictStructure` | `O(C + h·n)` (each candidate scored in context in O(n)) | `O(n + h + b)` | `C` = cost of `FindStemLoops`, `h` = number of candidate stem-loops, `b` = selected base pairs |
| `CalculateMfeStructure` / `PredictStructureMfe` (DP + traceback) | `O(n³)` fill + `O(n²)`-style traceback | `O(n²)` | Same recurrences as `CalculateMinimumFreeEnergy`; traceback re-evaluates options per visited cell |

A measured baseline is recorded in the MFE benchmark fixture ([RnaSecondaryStructure_MFE_Benchmark.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_MFE_Benchmark.cs)); the traceback adds an `O(n²)`-order pass over the already-filled `O(n³)` matrices and does not change the asymptotic class.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.CalculateMfeStructure(string, int)`: MFE-optimal structure via DP traceback; returns `MfeStructure` (sequence, dot-bracket, base-pair tuples, energy).
- `RnaSecondaryStructure.PredictStructureMfe(string, int)`: MFE-optimal structure in the `SecondaryStructure` shape (optimal counterpart of `PredictStructure`).
- `RnaSecondaryStructure.PredictStructure(string, int, int, int)`: stem-loop assembly (exact over the candidates, not the MFE fold).
- `RnaSecondaryStructure.FindStemLoops(string, int, int, int, bool)`: Supplies the candidate stem-loops used by `PredictStructure`.
- `RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int)`: Scalar global MFE routine documented in [RNA_Free_Energy.md](./RNA_Free_Energy.md); shares its `FillDp` matrices with the traceback.

### 5.2 Current Behavior

`CalculateMfeStructure` fills the V/W/WM matrices via the shared `FillDp` routine (identical to `CalculateMinimumFreeEnergy`, but into non-pooled arrays that survive the fill) and then runs the Zuker–Stiegler traceback to recover the optimal base-pair set. The structure energy is `W(0,n-1)` and equals the scalar MFE. `PredictStructureMfe` wraps this in the `SecondaryStructure` record (and, because the optimum is pseudoknot-free, its `Pseudoknots` field is empty).

`PredictStructure` scores each candidate stem-loop in context and runs the exterior recursion over the candidates; chosen spans are disjoint, so its `Pseudoknots` field (from `DetectPseudoknots(...)`) is always empty. Both paths emit dot-bracket strings using parentheses and dots only.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- RNA Watson-Crick and wobble pairing rules through the base-pair helper logic.
- Zuker–Stiegler (1981) MFE recurrences (V/W/WM) and the consistent traceback that recovers the optimal pseudoknot-free structure, with energy equal to the scalar MFE [2][6].
- Dot-bracket output for the selected secondary-structure pairing pattern.

**Intentionally simplified:**

- `PredictStructure` optimises only over combinations of isolated hairpin stem-loops; **consequence:** it is a restricted summary, retained alongside the optimal path. Use `CalculateMfeStructure` / `PredictStructureMfe` when global optimality is required.

**Not implemented:**

- Pseudoknotted optima (the standard `O(n³)` DP recurrences are pseudoknot-free by construction); **users should rely on:** dedicated pseudoknot folders for that structure class.
- Comparative or consensus secondary-structure prediction; **users should rely on:** no current alternative in this class.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input sequence | Returns empty strings, empty collections, and zero energy | Explicit early return in `PredictStructure` |
| No candidate stem-loops | Dot-bracket contains only `.` characters and `MinimumFreeEnergy = 0` | No selected stem-loops contribute base pairs or energy |
| Overlapping low-energy candidates | The combination of disjoint candidates with the lowest total in-context energy is chosen (ties: earlier candidate) | Exterior recursion over intervals |
| Candidate with ΔG ≥ 0 in context | Never selected (the open chain is better) | Zuker open-chain option W(j−1) |

### 6.2 Limitations

The MFE-optimal path returns the global optimum **for the Turner 2004 energy model and the pseudoknot-free structure class**: it does not predict pseudoknots (excluded by the `O(n³)` recurrences), modified bases, or comparative/consensus structure. The `PredictStructure` path is restricted to isolated hairpin stem-loops; its energy is the true energy of its structure but not the thermodynamic optimum; prefer the MFE path when optimality matters.

## 7. Examples and Related Material

- [RNA-STRUCT-001](../../../tests/TestSpecs/RNA-STRUCT-001.md) documents the repository's RNA secondary-structure test specification.
- [RnaSecondaryStructure_MfeStructure_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_MfeStructure_Tests.cs) is the canonical test fixture for the MFE-optimal traceback.
- [RNA_Stemloop.md](./RNA_Stemloop.md) documents the stem-loop finder used as the `PredictStructure` candidate generator.
- [RNA_Free_Energy.md](./RNA_Free_Energy.md) documents the underlying energy calculations and the scalar MFE.

## 8. References

1. Nussinov, R., and A. B. Jacobson. 1980. Fast algorithm for predicting the secondary structure of single-stranded RNA. Proceedings of the National Academy of Sciences 77(11):6309-6313.
2. Zuker, M., and P. Stiegler. 1981. Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. Nucleic Acids Research 9(1):133-148. https://pmc.ncbi.nlm.nih.gov/articles/PMC326673/
3. Mathews, D. H., et al. 2004. Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proceedings of the National Academy of Sciences 101(19):7287-7292.
4. Wikipedia contributors. Nucleic acid secondary structure. Wikipedia.
5. Wikipedia contributors. Nussinov algorithm. Wikipedia.
6. Washietl, S. (2012). 6.047/6.878 Lecture 08: RNA Folding. MIT. Fig. 13 (Zuker F/C/M/M¹ recurrences and traceback). https://web.mit.edu/6.047/book-2012/Lecture10_RNAStructure/Lecture10_RNAStructure_standalone.pdf

