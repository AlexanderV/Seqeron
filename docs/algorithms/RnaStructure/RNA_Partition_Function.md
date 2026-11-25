# RNA Partition Function (McCaskill, Turner 2004 model)

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-PARTITION-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-11-25 |

## 1. Overview

The McCaskill partition function computes the equilibrium statistical weight `Z` of the entire
ensemble of pseudoknot-free RNA secondary structures, and from it the equilibrium probability that
any two positions `(i,j)` form a base pair. Unlike minimum-free-energy folding, which returns a
single optimal structure, the partition function summarises the whole Boltzmann ensemble, so it
expresses folding uncertainty and structural alternatives [1]. `CalculatePartitionFunction` is the
counterpart of ViennaRNA 2.7.2 `fc.pf()` + `bpp()` + `mfe()` + `pr_structure()` +
`mean_bp_distance()` + `centroid()` on the **same Turner 2004 nearest-neighbour model** as the MFE
folder (d2 dangles); it replaces the earlier fixed-energy-per-pair teaching model (removed, see
§5.3). `CalculateStructureProbability` converts an energy gap into a Boltzmann probability.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In thermodynamic equilibrium each structure `S` occurs with probability
`p(S) = exp(−E(S)/RT) / Z`, `Z = Σ_S exp(−E(S)/RT)` [1][5]. Admissible pairs are Watson–Crick
(A–U, G–C) and the G·U wobble; hairpin loops enclose ≥ 3 nucleotides (NNDB).

### 2.2 Core Model

`β = 1/RT`, `RT = 1.98717·10⁻³ · T` kcal/mol (ViennaRNA `GASCONST`, `T` in K) [5].

- `Z = Σ_S exp(−βE(S))`; ensemble free energy `G = −RT·ln Z` (ViennaRNA `pf()` energy).
- `P(i,j) = (Σ_{S∋(i,j)} exp(−βE(S))) / Z`; `P_unpaired(i) = 1 − Σ_j P(i,j)`.
- `P(S) = exp((G − E(S))/RT)` (ViennaRNA `pr_structure`), in particular for the MFE structure.
- Ensemble diversity `⟨d⟩ = 2·Σ_{i<j} P(i,j)(1 − P(i,j))` (ViennaRNA `mean_bp_distance`); centroid =
  pairs with `P > 0.5` (minimal mean base-pair distance to the ensemble).

Inside recursion (unambiguous; every structure has exactly one derivation), with `b(E) = e^{−E/RT}`:

```
V(i,j)   = b(H(i,j)) + Σ_(p,q) b(I(i,j,p,q))·V(p,q) + GGUC corr. + b(Mclose(i,j))·Σ_k WM(i+1,k−1)·WM1(k,j−1)
WM1(i,j) = WM1(i,j−1) + V(i,j)·b(Mstem(i,j))                    (exactly one branch starting at i)
WM(i,j)  = Σ_k [ 1 + WM(i,k−1) ]·WM1(k,j)                        (≥ 1 branch, ML_base = 0)
W(j)     = W(j−1) + Σ_i W(i−1)·V(i,j)·b(Ext(i,j)),   Z = W(n)
```

`H, I, Mclose, Mstem, Ext` are the loop energies of the Turner 2004 model of `TurnerModel` (the
single evaluator shared with the MFE recursion, the traceback and `CalculateFreeEnergy`). Base-pair
probabilities come from McCaskill's **outside recursion** on the same grammar: `Ōv(i,j)` collects the
exterior loop `W(i−1)·b(Ext)·Wᵣ(j+1)`, every enclosing interior loop, the GGUC correction and the
multiloop branch / closing contributions (`Ōwm1`, `Ōwm`), and `P(i,j) = V(i,j)·Ōv(i,j)/Z` [1][3][5].

**Temperature.** `temperature` (K) rescales *every* Turner parameter with its enthalpy exactly as
ViennaRNA does (`params.c RESCALE_dG`): `ΔG(T) = ΔH − (ΔH − ΔG₃₇)·T/310.15`, with the ViennaRNA 2.7.2
Turner 2004 enthalpies (`Turner2004EnthalpyTables`); logarithmic loop extrapolations scale their
`1.75·R·T` / `1.08·T/T₃₇` prefactor; dangles and exterior/multiloop terminal mismatches are clamped to
≤ 0 (ViennaRNA: "must be <= 0"); NNDB-only terms with no published enthalpy (special G·U closure −2.2,
all-C hairpins, special-C bulge −0.9, GGUC/CUGG three-stack −4.12) stay at their 37 °C value. MFE,
traceback, structure evaluation and the partition function all read the same rescaled model, so
`G ≤ MFE` and `P(MFE)` hold at every temperature.

**Numerical stability.** Every inside/outside quantity of an `k`-nucleotide stretch is stored
multiplied by `us^k`, `us = exp(1.07·MFE/(RT·n))` (ViennaRNA `pf_scale` from the MFE, `sfact = 1.07`), so
long sequences neither overflow nor underflow; `Z` is returned as a double (+∞ beyond ≈ e^709) next to
its finite logarithm and `G`. When single-loop Boltzmann weights themselves leave the double range
(RT below ≈ 0.01 kcal/mol, a few Kelvin) the ensemble is the ground state (`G = MFE`, `P(MFE) = 1`).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Turner 2004 nearest-neighbour additivity, d2 dangles (ViennaRNA default). | Numeric `Z` inherits any model error; same model as the MFE folder. |
| ASM-02 | Pseudoknot-free structures only. | Crossing pairs are excluded from `Z`. |
| ASM-03 | Hairpin loops ≥ 3 nt (`minLoopSize` clamped to ≥ 3). | NNDB prohibits shorter loops. |
| ASM-04 | Rescaled parameters are used unrounded at `T ≠ 37 °C`; ViennaRNA's MFE truncates them to 0.01 kcal/mol (`int` cast), its partition function does not (`pf_smooth`). | `MFE` at `T ≠ 37 °C` differs from `fc.mfe()` by ≈ 0.05–0.1 kcal/mol (truncation noise), `G` by ≤ 0.03. |

### 2.4 Properties and Invariants

| ID | Invariant |
|----|-----------|
| INV-01 | `Z ≥ 1` (the open chain, `E = 0`, is one term); `G = −RT·ln Z ≤ MFE`. |
| INV-02 | Every `P(i,j) ∈ [0,1]`; rows of the probability matrix sum to ≤ 1. |
| INV-03 | `P_unpaired(i) + Σ_j P(i,j) = 1` at every position. |
| INV-04 | `P(MFE structure) = exp((G − E_mfe)/RT) ∈ (0,1]`. |
| INV-05 | `Z`, every `P(i,j)`, `MFE` and `P(MFE)` equal exhaustive enumeration over all structures scored by `CalculateFreeEnergy` (sequences ≤ 25 nt). |
| INV-06 | `CalculateStructureProbability(E, E) = 1`; `= exp((G−E)/RT)`; evaluated in the exponent (no overflow). |

### 2.5 Comparison with Related Methods

| Aspect | `CalculatePartitionFunction` | `CalculateMinimumFreeEnergy` / `CalculateMfeStructure` |
|--------|------------------------------|---------------------------------------------------------|
| Model | Turner 2004, d2, any `T` | the same |
| Output | `Z`, `P(i,j)`, `P_unpaired`, `G`, MFE + `P(MFE)`, diversity, centroid | single optimum |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | RNA (T read as U), case-insensitive; N / IUPAC / gaps never pair | not null; empty allowed |
| `minLoopSize` | `int` | 3 | minimum hairpin loop | values < 3 clamped to 3 |
| `temperature` | `double` | 310.15 | absolute temperature (K) | positive, finite |

### 3.2 Output / Return Value

`PartitionFunctionResult`: `PartitionFunction` (Z; +∞ beyond the double range), `BasePairProbabilities`
(0-based `(i<j)` → `P`, only `P > 0`), `EnsembleFreeEnergy` (`−RT·ln Z`, kcal/mol), `LogPartitionFunction`,
`UnpairedProbabilities`, `MfeStructure` / `MfeFreeEnergy` / `MfeProbability`, `EnsembleDiversity`,
`CentroidStructure`, `Temperature`. `CalculateStructureProbability(structureEnergy, ensembleEnergy,
temperature)` returns `exp((G − E)/RT)`. `CalculateUnpairedProbabilities` returns the legacy
`UnpairedProbabilityResult` projection of the same engine.

### 3.3 Preconditions and Validation

`null` → `ArgumentNullException`; `temperature ≤ 0`, NaN or ∞ → `ArgumentOutOfRangeException`; empty or
too-short input → `Z = 1`, all bases unpaired, no probabilities.

## 4. Algorithm

### 4.1 High-Level Steps

1. Fold the MFE (`FillDp` + traceback) at the requested temperature; derive the scale `us`.
2. Fill the scaled inside matrices `V, WM1, WM, W` (`FillPartitionInside`).
3. Fill the scaled outside matrices by decreasing span (`ComputePairProbabilities`); `P(i,j) = V·Ōv/Z`.
4. `G = −RT(ln Ẑ − n·ln us)`, unpaired, diversity, centroid, `P(MFE) = exp((G − E_mfe)/RT)`.

### 4.2 Decision Rules, Scoring, Reference Tables

Loop energies: `TurnerModel` (NNDB 37 °C tables `Turner2004Parameters`, `Turner2004InteriorLoopTables`;
enthalpies `Turner2004EnthalpyTables` generated from ViennaRNA 2.7.2 `default.c`, `intl11dH.h`,
`intl21dH.h`, `intl22dH.h`).

### 4.3 Complexity

`O(n³)` time (interior loops `O(n²·MAXLOOP²)`, MAXLOOP = 30), `O(n²)` memory (nine `n×n` double
arrays during the outside pass ≈ 290 MB at `n = 2000`). Measured (Release, one thread, MFE + inside + outside): `n = 200` 0.2 s,
`n = 500` 1.3 s, `n = 1000` 6.5 s, `n = 2000` 52 s (ViennaRNA 2.7.2: 0.9 / 4.6 / 29 s for 500 / 1000 / 2000).

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.PartitionFunction.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.PartitionFunction.cs),
[TurnerModel.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/TurnerModel.cs).

- `RnaSecondaryStructure.CalculatePartitionFunction(string, int minLoopSize, double temperature)`.
- `RnaSecondaryStructure.CalculateStructureProbability(double, double, double)`.
- `RnaSecondaryStructure.GenerateRandomRna(int, double)` / `(int, Random, double)`.
- MCP: none exposes the partition function (no wrapper to change).

### 5.2 Current Behavior

See §2.2. At 37 °C the model equals the NNDB table model bit for bit (loop terms rounded to 0.01 kcal/mol).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):** `Z`, `G`, `P(i,j)` (inside/outside), unpaired
probabilities, `P(S)`, mean base-pair distance and centroid [1][3][5]; temperature rescaling and `pf_scale` as ViennaRNA.

**Removed (was a simplification):** the fixed-energy-per-pair teaching model (Freiburg RNA tools) and its
`basePairEnergy` parameter: with the real Turner partition function implementable it only duplicated the engine.
Its combinatorial invariants (structure counts) are now verified by exhaustive enumeration on the Turner model.

**Not implemented:** pseudoknotted ensembles; ViennaRNA's `SMOOTH()` dangle weights (`pf_smooth`, ≤ 0.017 kcal/mol
per exterior end, specific-heat smoothing) — the exact Boltzmann sum over `eval_structure` energies is used instead.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status |
|---|------|------|--------|--------|
| 1 | ViennaRNA `pf()` weights exterior dangles with `SMOOTH()` | Reference quirk | `G` differs by ≤ 0.017 kcal/mol per sequence end | documented (Evidence §9) |
| 2 | Default ViennaRNA file rounds NNDB tables to 0.1 | Reference difference | `G` differs by ≈ 0.1 (mean) / 0.5 (max) from `fc.pf()` with the default file | documented |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior |
|------|-------------------|
| `null` | `ArgumentNullException` |
| empty / < `minLoopSize+2` nt / no admissible pair | `Z = 1`, `G = 0`, all unpaired, `P(MFE) = 1` |
| T, lower case | read as upper-case RNA |
| N, IUPAC, `-`, `.` | never pair |
| `temperature ≤ 0`, NaN, ∞ | `ArgumentOutOfRangeException` |
| very cold (< ≈ 5 K) | ground state: `G = MFE`, `P(MFE) = 1` |
| very long / cold | `Z = +∞`, `LogPartitionFunction` finite |

### 6.2 Limitations

Pseudoknot-free; nearest-neighbour model error is inherited; `O(n²)` memory limits practical length to a few
thousand nucleotides.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var r = RnaSecondaryStructure.CalculatePartitionFunction("GGGAAAACCC");
// r.MfeStructure "(((....)))", r.MfeFreeEnergy -2.02, r.EnsembleFreeEnergy -2.1535 (ViennaRNA pf -2.1529),
// r.MfeProbability 0.805, r.EnsembleDiversity 1.32, r.UnpairedProbabilities[0] 0.0592
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_PartitionFunction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_PartitionFunction_Tests.cs)
- Evidence: [RNA-PARTITION-001-Evidence.md](../../../docs/Evidence/RNA-PARTITION-001-Evidence.md)
- Related algorithms: [Turner_McCaskill_Partition_Function](./Turner_McCaskill_Partition_Function.md) (accessibility), [RNA_Free_Energy](./RNA_Free_Energy.md)

## 8. References

1. McCaskill, J. S. 1990. *Biopolymers* 29(6-7):1105-1119. https://doi.org/10.1002/bip.360290621 (PMID 1695107)
2. Will, S. 2011. McCaskill inside recursion, MIT 18.417. https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf
3. Will, S. 2011. McCaskill base-pair probabilities, MIT 18.417. https://math.mit.edu/classes/18.417/Slides/mccaskill2.pdf
4. Mathews, D. H. et al. 2004. *PNAS* 101:7287-7292 (Turner 2004); NNDB https://rna.urmc.rochester.edu/NNDB/turner04/
5. Lorenz, R. et al. 2011. ViennaRNA Package 2.0. *Algorithms Mol Biol* 6:26; ViennaRNA 2.7.2 sdist (`partfunc/`, `probabilities/`, `params/params.c`, `params/default.c`).
