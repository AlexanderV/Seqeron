# RNA Minimum Free Energy (MFE) Prediction

> **Baseline / reference method.** The classic Nussinov baseline documented here is a useful comparator; it uses a simplified energy model, not full thermodynamic fidelity (the Turner nearest-neighbor model is the higher-fidelity path). See [Legacy / Baseline Methods](../CANONICAL_MAP.md).

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-MFE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Computes the minimum free energy (ΔG°37, kcal/mol) of an RNA secondary structure by Zuker–Stiegler dynamic programming over the loop decomposition, scored with the Turner 2004 nearest-neighbor thermodynamic parameters [1][2]. The folding problem is recast as: over all pseudoknot-free secondary structures, find the one whose summed loop free energies is minimal. The result is exact for the implemented (affine-multiloop, no-coaxial-stacking) energy model — it is the global optimum of that model, not a heuristic. `CalculateMinimumFreeEnergy` returns the optimal score; `PredictStructure` returns a concrete base-pair set / dot-bracket using a greedy stem-loop assembly. Use it to estimate RNA stability and to identify favourable folds.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An RNA secondary structure is a set of non-crossing base pairs (A-U, G-C Watson-Crick, and G-U wobble). Anfinsen's thermodynamic hypothesis predicts the native fold is (near) the minimum-free-energy state. The nearest-neighbor model assigns free energy not to individual pairs but to the *loops* of the structure (stacked pairs, hairpin loops, bulge/interior loops, multibranch loops); the structure's ΔG° is the sum of loop contributions [1][2].

### 2.2 Core Model

Following Zuker & Stiegler (1981) [1], for subsequence indices `i < j` two dynamic-programming quantities are computed:

- **C(i,j)** — minimum free energy of the structure *closed* (enclosed) by pair (i,j). A closed structure is exactly one of: a hairpin loop, an interior/bulge loop over an inner pair, a stacked pair, or a multibranch loop [4]:

  `C(i,j) = min{ Hairpin(i,j), min_{i<i'<j'<j} [ Interior(i,j,i',j') + C(i',j') ], Multiloop(i,j) }`

- **W(j) / F** — minimum free energy of the unconstrained prefix/exterior segment, allowing the last base to be unpaired or to close a helix:

  `W(j) = min{ W(j−1), min_{0≤i≤j} [ W(i−1) + C(i,j) + end_penalty(i,j) + dangle(i,j) ] }`

- **WM / M, M1** — multibranch-loop fragment matrices accumulating ≥1 (M) and exactly-1 (M1) branch contributions [4].

The minimum hairpin loop is 3 nucleotides (`j − i − 1 ≥ 3`) [1][2]. Loop energies come from Turner 2004: nearest-neighbor stacking, hairpin/bulge/interior initiation tables, terminal AU/GU end penalties (+0.45), terminal-mismatch and first-mismatch bonuses, and special tri/tetra/hexaloop totals [2]. For a unimolecular fold the bimolecular helix-initiation constant is **not** applied [2, NNDB hairpin-example-1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Nearest-neighbor independence: a loop's energy depends only on its closing pair(s) and adjacent bases. | Non-local stabilization (tertiary contacts) is unmodeled; predicted ΔG° deviates from experiment. |
| ASM-02 | Fixed conditions: 37 °C, 1 M NaCl (Turner 2004 parameter conditions). | Different temperature/ionic strength shift all loop energies. |
| ASM-03 | Pseudoknot-free: only non-crossing pairs are scored. | Pseudoknotted folds are outside the optimization. |
| ASM-04 | Affine multiloop model with per-unpaired cost = 0 (offset 9.25, per-helix −0.63 from Turner 2004; unpaired coefficient simplified to 0). | Multiloop optima may shift relative to the full affine model; single-hairpin/stack/interior optima are unaffected. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | MFE ≤ 0 for every input. | The empty (open-chain) structure has ΔG° = 0 and is always in the search space; the minimum is therefore never positive [1]. |
| INV-02 | MFE is non-increasing under suffix extension: `MFE(prefix)` ≥ `MFE(prefix + suffix)`. | Extending the sequence only adds folding options; the optimum over a superset cannot be larger [1]. |
| INV-03 | The MFE score is deterministic: repeated evaluation of the same sequence yields the identical value. | The DP has no randomness; verified by the benchmark's determinism assertion. (The classic baseline uses a different, simplified energy model and is NOT expected to match the Turner-model score.) |

### 2.5 Comparison with Related Methods

| Aspect | MFE (this) | Partition function (McCaskill) |
|--------|------------|-------------------------------|
| Output | single optimal structure + ΔG° | base-pairing probabilities over the ensemble |
| Optimality | global optimum of the energy model | ensemble average |
| Complexity | O(n³) | O(n³) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `rnaSequence` | `string` | required | RNA sequence; case-insensitive, U-based | A/C/G/U; non-pairing chars never pair |
| `minLoopSize` | `int` | 3 | minimum hairpin loop length | clamped up to 3 if < 3 [1][2] |
| `minStemLength`/`maxLoopSize` | `int` | 3 / 10 | (`PredictStructure` only) greedy stem-loop search bounds | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (`CalculateMinimumFreeEnergy`) | `double` | optimal ΔG°37 in kcal/mol, rounded to 2 decimals (≤ 0) |
| (`PredictStructure`) `SecondaryStructure` | record | `Sequence`, `DotBracket`, `BasePairs`, `StemLoops`, `Pseudoknots`, `MinimumFreeEnergy` |

### 3.3 Preconditions and Validation

`null`/empty sequence → 0. Sequence shorter than `minLoopSize + 2` → 0 (cannot enclose a 3-nt hairpin). Input is upper-cased; indexing is 0-based; only A/C/G/U pair (Watson-Crick + G-U wobble); unrecognized characters simply never form pairs. No exceptions are thrown for unfoldable input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate/normalize input; short-circuit empty/too-short sequences to 0.
2. Fill `C(i,j)` and `WM(i,j)` by increasing span: evaluate hairpin, stacking, interior/bulge (bounded by MAXLOOP = 30), and multibranch options.
3. Accumulate the exterior optimum `W(j)` allowing unpaired bases, helix closure with AU/GU end penalties and dangling ends.
4. Return `W(n−1)` rounded to 2 decimals. `PredictStructure` instead greedily selects non-overlapping stem-loops by energy and emits base pairs + dot-bracket.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Turner 2004 parameter tables (kcal/mol, 37 °C) embedded in `RnaSecondaryStructure.cs` and sourced from NNDB [2]: nearest-neighbor stacking (WC + G-U), hairpin/bulge/interior loop initiation, terminal AU/GU end penalty (+0.45), terminal-mismatch stacking, first-mismatch bonuses (UU/GA −0.9, GG −0.8, special G-U closure −2.2), special tri/tetra/hexaloop totals, 1×1 interior loop table, dangling ends, and affine multibranch parameters (a = 9.25, c = −0.63). Internal/bulge loops are bounded by MAXLOOP = 30. A flat `double[]` buffer from `ArrayPool` backs the `n×n` matrices (O(n²) space).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateMinimumFreeEnergy` | O(n³) | O(n²) | interior loops bounded by MAXLOOP = 30 keep the interior search O(1) per (i,j) over a constant window; multiloop split is the O(n³) term [3][4]. |
| `PredictStructure` | O(n²) per stem-loop scan | O(n) | greedy heuristic, not the DP optimum. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int)`: Zuker-style O(n³) DP MFE score.
- `RnaSecondaryStructure.PredictStructure(string, int, int, int)`: greedy stem-loop structure + dot-bracket.
- `RnaSecondaryStructure.CalculateMinimumFreeEnergyClassic(string, int)`: internal O(n³) baseline using a simplified per-pair (Nussinov-style) energy model, retained only as a timing comparator in the benchmark — it does NOT use the Turner model and is not expected to match the Turner-model score.

### 5.2 Current Behavior

The DP buffers are rented from `ArrayPool<double>` and indexed as a flat `i*n + j` array. The interior/bulge search is bounded by `MAXLOOP = 30` and uses sliding lower bounds. A special GGUC/CUGG 3-stack context (NNDB note b) is detected and scored as a single −4.12 term. **Search reuse:** the suffix tree was evaluated and is **not** used — MFE is a thermodynamic-scoring dynamic program, not an exact-substring search; the suffix tree provides no benefit for energy minimization (recorded per "Reuse existing infrastructure"). `PredictStructure` is a greedy heuristic and may return a structure whose energy is above the true MFE; for the optimal score use `CalculateMinimumFreeEnergy`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Zuker–Stiegler loop decomposition (hairpin / stacking / bulge-interior / multibranch) with C(i,j) and exterior W recurrences [1][4].
- Turner 2004 nearest-neighbor stacking, hairpin/bulge/interior initiation, terminal AU/GU end penalty (+0.45), terminal-mismatch and first-mismatch bonuses, special tri/tetra/hexaloop totals [2].
- Minimum hairpin loop = 3 nt; no helix-initiation constant for unimolecular folds [1][2].

**Intentionally simplified:**

- Multibranch loop: affine model with per-unpaired-nucleotide cost set to 0 (ASM-04); **consequence:** multiloop optima can differ from the full affine Turner model, though single-hairpin/stack/interior optima (the validated worked examples) are unaffected.
- Coaxial stacking is not applied inside multiloops/exterior junctions; **consequence:** energies of multi-helix junctions are slightly less stable than ViennaRNA/RNAstructure defaults.

**Not implemented:**

- Pseudoknots; **users should rely on:** `DetectPseudoknots` for detection only — pseudoknotted MFE is out of scope (ASM-03).
- Logarithmic multiloop model and ensemble/partition-function quantities; **users should rely on:** an external tool (ViennaRNA/RNAstructure) for those.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Multiloop per-unpaired cost = 0 | Assumption | shifts multiloop optima only | accepted | ASM-04; single-hairpin path validated against NNDB examples |
| 2 | 2-decimal rounding vs NNDB 1-decimal | Assumption | display precision | accepted | tests assert `.Within(1e-2)` vs NNDB total and `.Within(1e-9)` vs exact sum |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` / `""` | 0 | no sequence |
| Homopolymer (`AAAAAAAA`) | 0 | no base pair possible; open chain ΔG = 0 [1] |
| Length < `minLoopSize + 2` | 0 | cannot enclose a 3-nt hairpin [1][2] |
| `minLoopSize < 3` | clamped to 3 | nearest-neighbor rules prohibit shorter loops [2] |

### 6.2 Limitations

Pseudoknot-free only; no coaxial stacking; simplified multiloop unpaired cost; fixed 37 °C / 1 M NaCl Turner conditions; `PredictStructure` is a greedy heuristic and not guaranteed to reach the DP optimum. Not validated against a published numeric multiloop optimum (single-hairpin worked examples are the cited evidence).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// NNDB Turner 2004 Hairpin Example 1 reconstructed as a full hairpin:
// stem CACA / UGUG closing pairs C-G,A-U,C-G,A-U; 6-nt loop AAAAAA.
double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy("CACAAAAAAAUGUG");
// mfe == -1.41  (NNDB tabulated -1.4)
var s = RnaSecondaryStructure.PredictStructure("CACAAAAAAAUGUG");
// s.DotBracket == "((((......))))", s.BasePairs.Count == 4
```

**Numerical walk-through (NNDB Example 1 [5]):**
3 stacks (CG/AU −2.11, AU/CG −2.24, CG/AU −2.11) + AU end penalty (+0.45) + terminal mismatch AU·AA (−0.8) + hairpin initiation(6) (+5.4) = **−1.41 kcal/mol**.

### 7.2 Performance Baseline

`CalculateMinimumFreeEnergy` (full Turner model), measured via `RnaSecondaryStructure_MFE_Benchmark` (median of 3 runs, Release, Apple-silicon dev machine, 2026-06-14):

| Length n | Time (ms) | Notes |
|----------|-----------|-------|
| 100 | ~25 | — |
| 200 | ~83 | ×3.3 for ×2 length |
| 300 | ~214 | — |
| 500 | ~639 | — |
| 1000 | ~3155 | ×5 for ×2 length |

The ~cubic growth (≈ ×8 per length-doubling for the dominant term, partially offset by the bounded MAXLOOP = 30 interior search) is consistent with the O(n³) complexity in §4.3.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_MinimumFreeEnergy_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_MinimumFreeEnergy_Tests.cs) — covers `INV-01`, `INV-02`
- Performance baseline: [RnaSecondaryStructure_MFE_Benchmark.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_MFE_Benchmark.cs) — `INV-03` + timing
- Evidence: [RNA-MFE-001-Evidence.md](../../../docs/Evidence/RNA-MFE-001-Evidence.md)
- Related algorithms: [Hairpin_Energy_Calculation](Hairpin_Energy_Calculation.md), [RNA_Base_Pairing](RNA_Base_Pairing.md)

### 7.4 Change History

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-14 | 1.0 | Initial documentation (RNA-MFE-001). |

## 8. References

1. Zuker M, Stiegler P. 1981. Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. Nucleic Acids Research 9(1):133–148. https://doi.org/10.1093/nar/9.1.133
2. Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH. 2004. Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. PNAS 101:7287–7292. https://doi.org/10.1073/pnas.0401799101
3. Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. 2011. ViennaRNA Package 2.0. Algorithms for Molecular Biology 6:26. https://doi.org/10.1186/1748-7188-6-26
4. Ward M, Datta A, Wise M, Mathews DH. 2017. Advanced multi-loop algorithms for RNA secondary structure prediction reveal that the simplest model is best. Nucleic Acids Research 45(14):8541–8552. https://doi.org/10.1093/nar/gkx512
5. NNDB Turner 2004 Hairpin Loop Example 1. https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-1.html
