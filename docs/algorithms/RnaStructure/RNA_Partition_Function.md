# RNA Partition Function (McCaskill)

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-PARTITION-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

The McCaskill partition function computes the equilibrium statistical weight `Z` of the
entire ensemble of pseudoknot-free RNA secondary structures, and from it the equilibrium
probability that any two positions `(i,j)` form a base pair. Unlike minimum-free-energy
folding, which returns a single optimal structure, the partition function summarises the
whole Boltzmann ensemble, so it can express folding uncertainty and structural
alternatives [1]. The computation is an exact dynamic program of order `O(n³)` time and
`O(n²)` space [1][2]. This implementation uses a *simplified* fixed-per-base-pair energy
model (the model of the Freiburg RNA teaching tool) [4]; the partition-function
mathematics and all derived probabilities are exact for that energy model.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An RNA secondary structure is a set of non-crossing (pseudoknot-free) base pairs over a
sequence. In thermodynamic equilibrium each structure `S` occurs with Boltzmann
probability `p(S) = exp(−E(S)/RT) / Z`, where `Z` is the normalising partition function
[1][5]. Admissible pairs are Watson-Crick (A-U, G-C) and the GU wobble pair [4].

### 2.2 Core Model

Let `β = 1/RT`. The partition function over all structures is [1][4]:

```
Z = Σ_S exp(−E(S)/RT)
```

McCaskill's inside dynamic program uses two matrices over 0-based inclusive
sub-sequences `[i..j]` [2]:

- `Q[i,j]` = partition function of `[i..j]`; the empty interval (`i > j`) is defined as 1.
- `Qᵇ[i,j]` = partition function of `[i..j]` restricted to structures in which `(i,j)` is paired.

Recurrences (with minimum hairpin loop `m`) [2]:

```
Q[i,j]  = Q[i,j-1] + Σ_{i ≤ k < j-m} Q[i,k-1] · Qᵇ[k,j]          (Q[i,j] = 1 for i ≥ j-m)
Qᵇ[i,j] = exp(−β·E_bp) · Q[i+1,j-1]   if (i,j) can pair and j-i > m, else 0
```

The total partition function is `Z = Q[1,n]` [2]. The equilibrium base-pair probability is
the **outside recursion** [3]: `P[i,j] = Qᵇ[i,j]·O[i,j]/Z`, where the outside partition
function `O[i,j]` collects every structure of the sequence outside `[i..j]` together with
every enclosing pair:

```
O[i,j] = Q[1,i-1]·Q[j+1,n]                                              (external: (i,j) unenclosed)
       + Σ_{k<i, l>j, CanPair(k,l), l-k>m}  w·Q[k+1,i-1]·Q[j+1,l-1]·O[k,l]   (enclosing pairs)
P[i,j] = Qᵇ[i,j] · O[i,j] / Z
```

with `w = exp(−β·E_bp)`. The **external term alone is not sufficient**, even in the flat
fixed-per-pair model: a pair that can be nested inside another pair receives a strictly larger
probability from the enclosing terms (e.g. `GGGAAACCC`, E_bp=0: P[2,6]=6/20=0.30, whereas the
external term gives only 1/20=0.05).

The single-structure Boltzmann probability is `p(S) = exp(−βE(S)) / Z` [1][5].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Each base pair contributes a single fixed energy `E_bp`; loop-type (hairpin/internal/bulge/multiloop), stacking and dangling-end energies are ignored | Numeric `Z` and per-pair probabilities differ from a full Turner 2004 model (e.g. ViennaRNA); rankings/invariants are unaffected |
| ASM-02 | Structures are pseudoknot-free (non-crossing pairs only) | Crossing/pseudoknotted ensembles are not counted [1] |
| ASM-03 | Minimum hairpin loop = 3 unpaired bases | Pairs with span ≤ 3 are excluded (steric constraint) [stdmodel] |
| ASM-04 | Fixed temperature (default 310.15 K) and `R = 1.987 cal/(mol·K)` | `Z` scales with `RT`; default matches the Turner/NNDB 37 °C convention [5] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Z ≥ 1` | the empty structure always contributes `exp(0) = 1` (base case `Q = 1`) [2] |
| INV-02 | `P[i,j] ∈ [0,1]` | `P[i,j]` is a sum of Boltzmann probabilities of structures containing `(i,j)` [1][3] |
| INV-03 | If `E_bp = 0`, `Z` = number of admissible structures | every Boltzmann weight is `exp(0) = 1` [4] |
| INV-04 | `Z` strictly increases as `E_bp` decreases | each pair weight `exp(−β·E_bp)` increases, and `Z` is a positive combination of them [1] |
| INV-05 | `p(S) = exp(−βE(S))/Z ∈ (0,1]`, `= 1` when `S` is the whole ensemble | normalised Boltzmann factor [1][5] |

### 2.5 Comparison with Related Methods

| Aspect | Partition function (McCaskill) | MFE folding (Zuker) |
|--------|-------------------------------|---------------------|
| Output | ensemble weight `Z` + all base-pair probabilities | one optimal structure + its ΔG |
| Captures alternatives | yes (full Boltzmann ensemble) | no |
| Complexity | O(n³) time, O(n²) space | O(n³) time, O(n²) space |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | RNA sequence | A/C/G/U, case-insensitive; T treated as U; null → exception |
| `basePairEnergy` | `double` | `-1.0` kcal/mol | fixed `E_bp` per base pair | any real value |
| `temperature` | `double` | `310.15` K | absolute temperature | must be > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `PartitionFunction` | `double` | `Z = Σ_S exp(−E(S)/RT)`; `Z ≥ 1` |
| `BasePairProbabilities` | `IReadOnlyDictionary<(int I,int J),double>` | `P[i,j]` for every pair (0-based, `i<j`) that can form |

`CalculateStructureProbability(structureEnergy, ensembleEnergy, temperature)` returns the
scalar Boltzmann probability `exp(−βE_S)/exp(−βE_ens)`.

### 3.3 Preconditions and Validation

Null `sequence` → `ArgumentNullException`. Non-positive `temperature` →
`ArgumentOutOfRangeException`. Empty `sequence` → `Z = 1`, empty probability map
(only the empty structure). Indexing is 0-based; base-pair keys satisfy `i < j`.
The sequence is upper-cased internally; only A-U, G-C, and GU pairs are admissible.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; handle empty sequence (`Z = 1`).
2. Compute `RT` and the per-pair Boltzmann weight `exp(−E_bp/RT)`.
3. Fill `Qᵇ` and `Q` by increasing sub-sequence length (inside recursion).
4. Read `Z = Q[0,n-1]`.
5. Compute every base-pair probability via the outside recursion
   `P[i,j] = Qᵇ[i,j]·O[i,j]/Z` (external term plus enclosing-pair terms; see §2.2),
   processing pairs outermost-first so each enclosing `O[k,l]` is ready.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Pairing table: A-U, G-C (Watson-Crick) and G-U (wobble) only, via `CanPair` [4].
- Minimum hairpin loop `m = 3`: a pair `(i,j)` requires `j − i > 3`.
- Per-pair energy `E_bp` (default −1.0 kcal/mol); `R = 1.987 cal/(mol·K)`; default `T = 310.15 K` [5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculatePartitionFunction` | O(n³) | O(n²) | two `n×n` matrices; triple loop over `i,j,k` [1][2] |
| `CalculateStructureProbability` | O(1) | O(1) | two exponentials |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.CalculatePartitionFunction(string, double, double)`: McCaskill `Z` + base-pair probabilities.
- `RnaSecondaryStructure.CalculateStructureProbability(double, double, double)`: Boltzmann probability `exp(−βE)/Z`.
- `RnaSecondaryStructure.GenerateRandomRna(int, double)` / `(int, Random, double)`: random RNA generation (seeded overload is deterministic).

### 5.2 Current Behavior

The DP fills matrices by increasing interval length so all dependencies precede each cell.
Base-pair probabilities are computed via the full **outside recursion** `P[i,j] = Qᵇ[i,j]·O[i,j]/Z`
(see §2.2): the external term plus, for every enclosing pair, `w·Q[k+1,i-1]·Q[j+1,l-1]·O[k,l]`.
Pairs are processed outermost-first so each enclosing `O[k,l]` is already known. (An earlier
version computed only the external term; that under-reported the probability of any nestable
pair and was corrected on 2026-06-16 — see the validation report.) The implementation does
**not** use the repository suffix tree: this is a numerical
dynamic program over base-pairing weights, not a substring-search task, so the suffix tree
does not apply (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `Z = Σ_S exp(−E(S)/RT)` and `Z = Q[1,n]` [1][2].
- Inside recursion `Q[i,j] = Q[i,j-1] + Σ_k Q[i,k-1]·Qᵇ[k,j]`, base case `Q = 1` [2].
- `Qᵇ[i,j] = exp(−β E_bp)·Q[i+1,j-1]` for admissible pairs [2][4].
- Base-pair probability via the outside recursion `P[i,j] = Qᵇ[i,j]·O[i,j]/Z`, `O` per §2.2 [3].
- Boltzmann structure probability `p = exp(−βE)/Z` [1][5].

**Intentionally simplified:**

- Energy model: fixed per-pair `E_bp` instead of Turner 2004 loop/stacking/dangling energies (ASM-01); **consequence:** numeric `Z` and probabilities differ from ViennaRNA on real sequences, though invariants and ensemble interpretation are preserved.

**Not implemented:**

- Pseudoknotted ensembles; **users should rely on:** no current alternative in-repo (out of scope, ASM-02).
- Full Turner-parameter partition function (interior/multiloop Boltzmann terms); **users should rely on:** ViennaRNA `RNAfold -p` for production thermodynamic ensembles.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` sequence | `ArgumentNullException` | contract |
| empty sequence | `Z = 1`, no pairs | only the empty structure exists [2] |
| no admissible pair (`"AAAA"`) | `Z = 1` | INV-01 |
| pair span ≤ 3 only (`"GC"`) | `Z = 1` | min hairpin loop (ASM-03) |
| `temperature ≤ 0` | `ArgumentOutOfRangeException` | physically invalid |

### 6.2 Limitations

Only Watson-Crick + GU pairing; pseudoknot-free; simplified energy model (ASM-01); not
intended as a substitute for a full Turner-parameter partition function. For large `n` the
`O(n³)` time and `O(n²)` memory dominate.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var result = RnaSecondaryStructure.CalculatePartitionFunction("GGGAAACCC", basePairEnergy: 0.0);
// result.PartitionFunction == 20.0  (number of admissible structures, since E_bp = 0)
// result.BasePairProbabilities[(0,8)] == 0.30
```

**Numerical walk-through:** with `E_bp = 0`, every Boltzmann weight is 1, so
`Z("GGGAAACCC") = 20` is exactly the number of non-crossing base-disjoint structures.
6 of those 20 structures contain the pair (0,8), hence `P[0,8] = 6/20 = 0.30`.

**Performance baseline (Phase 8):** `CalculatePartitionFunction` on a random length-300
RNA (single-threaded, .NET 9, Apple Silicon, Release) completes well under 1 second;
the O(n³) growth is the dominant factor. Measured during validation — see §7.3 tests for
the property/performance test.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_PartitionFunction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PartitionFunction_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [RNA-PARTITION-001-Evidence.md](../../../docs/Evidence/RNA-PARTITION-001-Evidence.md)
- Related algorithms: [RNA_Free_Energy](./RNA_Free_Energy.md)

## 8. References

1. McCaskill, J. S. 1990. The equilibrium partition function and base pair binding probabilities for RNA secondary structure. *Biopolymers* 29(6-7):1105-1119. https://doi.org/10.1002/bip.360290621 (PMID 1695107)
2. Will, S. 2011. McCaskill algorithm (inside recursion) lecture notes, MIT 18.417. https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf
3. Will, S. 2011. McCaskill base-pair probabilities (outside recursion) lecture notes, MIT 18.417. https://math.mit.edu/classes/18.417/Slides/mccaskill2.pdf
4. Freiburg RNA Tools. McCaskill teaching tool (simplified fixed-per-pair energy model). https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill
5. ViennaRNA Package. Partition Function and Equilibrium Properties (pf_fold reference). https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/pf_fold.html
