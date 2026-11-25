# Turner-2004 McCaskill Partition Function and Structural Accessibility (SA)

| Field | Value |
|-------|-------|
| Algorithm Group | RNA_Structure |
| Test Unit ID | MIRNA-TARGET-001 (SA wiring); new RNA-STRUCT capability |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Annotation |
| Implementation Status | Production |
| Last Reviewed | 2026-06-25 |

## 1. Overview

This is the McCaskill (1990) equilibrium **partition function** evaluated on the SAME Turner-2004
nearest-neighbour energy model used by the library's Zuker–Stiegler MFE folder. Instead of taking
the minimum free energy over all pseudoknot-free secondary structures, it sums the Boltzmann
weight `exp(−E(S)/RT)` over the whole ensemble, yielding the partition function `Z`, the
equilibrium base-pair probabilities `P(i,j)`, the per-base unpaired probabilities
`p_unpaired(i) = 1 − Σ_j P(i,j)`, and the ensemble free energy `−RT·ln Z`. It also computes the
**joint region-unpaired probability** (RNAplfold-style accessibility) used by the TargetScan
context++ **SA** feature. The partition function is exact for the Turner-2004 model; the result
is deterministic. [1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An RNA molecule does not adopt a single structure but a Boltzmann-weighted ensemble. McCaskill
showed the equilibrium partition function and base-pair binding probabilities can be computed by a
dynamic program with the same `O(N³)`/`O(N²)` complexity as MFE folding. [1] miRNA target-site
efficacy depends on how *accessible* (unpaired) the 3'UTR site is, which TargetScan scores from
RNAplfold unpaired probabilities (the SA feature). [3][4]

### 2.2 Core Model

For a sequence `S`, over the set `Ω` of pseudoknot-free structures with energy `E(s)` under the
Turner-2004 model: [1][2]

- Partition function: `Z = Σ_{s∈Ω} exp(−β E(s))`, with `β = 1/(RT)`, `RT = R·T`, `R = 1.98717×10⁻³ kcal/(mol·K)` (ViennaRNA `GASCONST`).
- Boltzmann probability of a structure: `p(s) = exp(−β E(s)) / Z`.
- Base-pair probability: `P(i,j) = (Σ_{s ∋ (i,j)} exp(−β E(s))) / Z` (inside/outside recursion). [1]
- Per-base unpaired probability: `p_unpaired(i) = 1 − Σ_j P(i,j)`. [2]
- Ensemble free energy: `ΔG_ensemble = −RT·ln Z`.
- Region accessibility: `P([a..b] unpaired) = Z_open(a,b) / Z`, where `Z_open(a,b)` sums the
  weights of structures in which no base of `[a..b]` is paired. RNAplfold's `_lunp` row `i`,
  column `L`, reports `P([i−L+1..i] unpaired)`. [3][4]

The inside recursion mirrors the MFE Zuker recurrence with `min → Σ` and `+ΔG → × exp(−ΔG/RT)`:
`Vexp(i,j)` (i·j paired) sums hairpin / stacking / internal / bulge / multibranch closings using
the SAME Turner loop-energy functions; `WMexp` is the multibranch region; `Wexp` is the external
loop. `Z = Wexp(0,n−1)`. [1]

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Nearest-neighbour additivity of the Turner-2004 model (loop energies independent). | Same modelling assumption as the MFE folder; ensemble probabilities inherit any model error. |
| ASM-02 | Pseudoknot-free structures only. | Crossing (pseudoknotted) pairs are excluded from `Z` and from accessibility. |
| ASM-03 | Temperature `T` (default 310.15 K): every Turner-2004 parameter is rescaled with its enthalpy, `ΔG(T) = ΔH − (ΔH − ΔG₃₇)·T/310.15` (ViennaRNA), not only `RT`. | Probabilities are temperature-specific; see [RNA_Partition_Function](./RNA_Partition_Function.md) §2.2. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Z ≥ 1` (`> 0` always). | The open chain (`E = 0`) is always one term, contributing `exp(0) = 1`. |
| INV-02 | `P(i,j) ∈ [0,1]` for every pair. | A ratio of a non-negative sub-sum to the total `Z`. |
| INV-03 | `Σ_j P(i,j) ≤ 1` for each `i`. | At most one pair can involve position `i` in any nested structure. |
| INV-04 | `p_unpaired(i) = 1 − Σ_j P(i,j) ∈ [0,1]`. | Complement of the paired probability. [2] |
| INV-05 | `ΔG_ensemble = −RT·ln Z ≤ MFE`. | The MFE structure is one term in `Z`, so `Z ≥ exp(−MFE/RT)`. [2] |
| INV-06 | Region accessibility `Z_open/Z ∈ [0,1]`. | `Z_open ≤ Z` (its structures are a subset). |

### 2.5 Comparison with Related Methods

The full ensemble result (`CalculatePartitionFunction`: `Z`, `P(i,j)`, `p_unpaired`, `ΔG_ens`, MFE structure and probability,
diversity, centroid) is documented in [RNA_Partition_Function](./RNA_Partition_Function.md); `CalculateUnpairedProbabilities`
returns the same engine's `Z`, `P(i,j)`, `p_unpaired`, `ΔG_ens`. The MFE folder (`CalculateMinimumFreeEnergy`) is the single-structure counterpart.

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `rnaSequence` | `string` | required | RNA (or DNA; T→U) sequence | case-insensitive; null/empty allowed |
| `minLoopSize` | `int` | 3 | minimum hairpin loop | clamped to ≥3 (NNDB) |
| `temperature` | `double` | 310.15 | absolute temperature (K) | positive, finite |
| `windowEnd` / `windowLength` | `int` | — | accessibility window (0-based inclusive end, length u) | window must fit in the sequence; length > 0 (signature `CalculateRegionUnpairedProbability(seq, windowEnd, windowLength, minLoopSize, temperature)`) |

### 3.2 Output / Return Value

`UnpairedProbabilityResult`:

| Field | Type | Description |
|-------|------|-------------|
| `PartitionFunction` | `double` | `Z` |
| `BasePairProbabilities` | `IReadOnlyDictionary<(int,int),double>` | `P(i,j)` for 0-based `i<j` (only `P>0`) |
| `UnpairedProbabilities` | `IReadOnlyList<double>` | `p_unpaired(i)`, one per position |
| `EnsembleFreeEnergy` | `double` | `−RT·ln Z` (kcal/mol) |

`CalculateRegionUnpairedProbability` returns `double` — `P(window unpaired) ∈ [0,1]`; `CalculateRegionOpeningEnergy` returns the opening energy `−RT·ln P` (kcal/mol; 0 when `P = 1`, `+∞` when `P` underflows to 0; Kertesz et al. 2007). MCP: `rna_region_accessibility`.

### 3.3 Preconditions and Validation

0-based inclusive coordinates. DNA accepted (T read as U); upper-cased internally. Null/empty →
`Z = 1`, no pairs, `ΔG = 0`. Sequence too short to pair → `Z = 1`, all bases unpaired. A
non-positive temperature throws `ArgumentOutOfRangeException`. An accessibility window that does
not fit throws `ArgumentOutOfRangeException` (this includes a null/empty sequence, a negative or ≥ n `windowEnd`, a window longer than the sequence and `windowLength ≤ 0`; unlike the per-base method, the region method has no "unit partition" result for an empty sequence); a window in a sequence too short to pair returns 1.0. `N`/IUPAC bases never pair, so a region of such bases has probability exactly 1. A `minLoopSize` > 3 removes the short hairpins from both `Z` and `Z_open`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Fold the MFE and fill the scaled inside matrices (`FillPartitionInside`), `Z = W(n)`; fill the outside matrices
   (`ComputePairProbabilities`) → `P(i,j)`; `p_unpaired(i) = 1 − Σ_j P(i,j)` (see [RNA_Partition_Function](./RNA_Partition_Function.md)).
2. `ΔG_ensemble = −RT·ln Z`.
3. Accessibility: `Z_open(window)/Z` — the same inside fill with every window base forbidden from pairing (both fills share one scale,
   which cancels in the ratio; the scale is ViennaRNA's default estimate `−185 + (T−37)·7.27` cal/(mol·nt) since no MFE is needed).

### 4.2 Decision Rules, Scoring, Reference Tables

The energy contributions are the Turner-2004 NNDB tables already used by the MFE folder
(stacking, hairpin/internal/bulge initiation and mismatches, int11/int21/int22, d2 terminal
mismatches/dangles, terminal-AU, linear multibranch `a/b/c` = 9.3/0/−0.9). [5] The inside recursion
is unambiguous (single-branch matrix WM1 as in ViennaRNA's fM1; each structure counted once); with
NNDB-precision parameters the ensemble free energy equals ViennaRNA 2.7.2 `fc.pf()` (d2) to
≤ 0.017 kcal/mol (mean 0.0004; 300 random sequences 10–150 nt, no repeated neighbours) — [Evidence RNA-PARTITION-001 §9](../../Evidence/RNA-PARTITION-001-Evidence.md). Each is converted to a Boltzmann factor `exp(−ΔG/RT)`. The **SA** feature
(`MiRnaAnalyzer`) reads, for a seed-matched site, the accessibility of the 14-nt window ending 7
nt 3' of the seed-match start (7mer-A1 decrements the start by 1), `log10`-transforms it, and
min-max scales by the verbatim Agarwal-2015 SA coefficients. RNAplfold reference parameters:
`-L 40 -W 80 -u 20`. [3][4]

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Z`, `p_unpaired`, `P(i,j)` | `O(n³)` | `O(n²)` | inside + outside recursion (McCaskill), no per-pair refold |
| `CalculateRegionUnpairedProbability` | `O(n³)` | `O(n²)` | two inside fills (full + window-forbidden) |
| SA over a 3'UTR site | `O(W³)` | `O(W²)` | local fold of a `W=80` context |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.PartitionFunction.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.PartitionFunction.cs)

- `RnaSecondaryStructure.CalculateUnpairedProbabilities(seq, minLoopSize, temperature)`: `Z`, `P(i,j)`, `p_unpaired`, `ΔG_ensemble` (projection of `CalculatePartitionFunction`).
- `RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd, windowLength, …)`: joint region-unpaired probability (`Z_open/Z`).
- `MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(...)` → `SaContribution`: SA feature ([MiRnaAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)).

### 5.2 Current Behavior

`P(i,j)` comes from McCaskill's outside recursion, `p_unpaired(i) = 1 − Σ_j P(i,j)`, so `Σ_j P(i,j) + p_unpaired(i) = 1` to
floating-point precision (exhaustive-enumeration test: 1e-9). SA folds a local `W=80` context around the window (RNAplfold's
local-folding intent; a base pairs only within ±`L=40`). When the SA window does not fit the UTR, SA contributes 0 and is reported
in `OmittedFeatures` (matching the perl's missing-`_lunp` → 0).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `Z = Σ exp(−E/RT)` on the Turner-2004 model; `P(i,j)`, `p_unpaired(i) = 1 − Σ_j P(i,j)`, `ΔG_ensemble = −RT·ln Z`. [1][2]
- Region accessibility `Z_open/Z` and the TargetScan SA window (row `utrStart+7`, `L=14`, `log10`, min-max scaled by the SA row). [3][4]

**Scope of the accessibility value (not a simplification of the unit):** `CalculateRegionUnpairedProbability` is the WHOLE-sequence ensemble quantity `Z_open/Z`, i.e. what RNAplfold `-u` reports for window `W ≥ n` and span `L ≥ n`. Verified against ViennaRNA 2.7.2: exact hard-constraint ratio `exp(−(G_constrained − G_ens)/RT)` (1440 random regions u = 1…n, 8–200 nt, 20/37/60 °C: max |Δ| 0.0056, mean 2e-5; the maximum is ViennaRNA's `SMOOTH()` end-dangle quirk) and `RNA.pfl_fold_up(seq, u, n, n)` (600 random queries: max |Δ| 1.0e-3; identical to the constraint ratio to 2e-6 on the locked cases). The joint probability is **not** the product of per-base `p_unpaired` (e.g. `GGGGAAAACCCCUUUUGGGG`, 5′ arm: joint > 1.2 × product).

**Intentionally simplified:**

- Local SA fold uses a single `W=80` context window rather than RNAplfold's full sliding-window average over all `W`-windows; **consequence:** the accessibility of the 14-nt site is captured (bases pair only within ±`L=40`), but values can differ slightly from a full RNAplfold run near long-range context.

**Not implemented:**

- RNAplfold's sliding-window local accessibility (`W < n`, maximum pair span `L < n`, Bernhart et al. 2006; ViennaRNA `pf_window.c`, `compute_pU`). It is not the arithmetic mean of independent per-window ensembles: the mean of `Z_open/Z` over the windows `[a, a+W−1]` that contain the region, each folded with `md.max_bp_span = L`, differs from `RNA.pfl_fold_up(seq, u, W, L)` by up to 0.07 (60 nt, W = 30, L = 15), so it needs the exterior-context bookkeeping of `pf_window.c` (≈ 1500 lines) to be ported. Registered as a separate work item; until then the exact RNAplfold value can be supplied to `ScoreTargetSiteContextPlusPlus` (`PlfoldUnpairedProbability`).
- Pseudoknotted contributions to `Z` / accessibility; **users should rely on:** no nearest-neighbour ensemble model includes them (energy-model floor).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Single `W=80` SA window instead of RNAplfold's sliding-window average | Simplification | see §5.3 | accepted | MiRNA unit (B13) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty | `Z = 1`, no pairs, `ΔG = 0` | only the empty structure |
| too short to pair | `Z = 1`, all bases unpaired | no admissible pair |
| `temperature ≤ 0` | `ArgumentOutOfRangeException` | `RT` undefined / non-physical |
| SA window off the UTR end | SA = 0, listed in `OmittedFeatures` | matches missing-`_lunp` → 0 |

### 6.2 Limitations

Pseudoknot-free only; nearest-neighbour model error is inherited; the SA local-fold is a single
`W=80` window, not the full RNAplfold sliding-window average. PCT (conservation) is not part of
this capability.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (the analytic tiny case).** `GAAAC` can adopt only the open chain or the
single hairpin closed by G(0)–C(4) over the 3-nt loop `AAA` (Turner hairpin energy 5.4 kcal/mol,
G–C ⇒ no terminal-AU penalty). With `RT = 1.98717·310.15/1000 = 0.6163207755`:

- `Z = 1 + exp(−5.4/RT) = 1.0001566226391072`
- `P(0,4) = exp(−5.4/RT)/Z = 0.00015659811229750914`
- `p_unpaired(0) = p_unpaired(4) = 0.9998434018877025`; `p_unpaired(2) = 1`
- `ΔG_ensemble = −RT·ln Z ≈ −9.652×10⁻⁵ ≤ MFE = 0`
- `CalculateRegionUnpairedProbability("GAAAC", 4, 5) = 1/Z` (only the open chain leaves the whole window unpaired).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_UnpairedProbabilities_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_UnpairedProbabilities_Tests.cs) — covers `INV-01`..`INV-06` and the analytic case.
- Tests: [MiRnaAnalyzer_TargetPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/MiRnaAnalyzer_TargetPrediction_Tests.cs) — `CTX-SA-001` (SA wiring).
- Evidence: [MIRNA-TARGET-001-Evidence.md](../../../docs/Evidence/MIRNA-TARGET-001-Evidence.md) (§ SA — structural accessibility).
- Related algorithms: [RNA_Partition_Function](./RNA_Partition_Function.md) (partition-function unit), [RNA_Free_Energy](./RNA_Free_Energy.md) (MFE folder).

## 8. References

1. McCaskill, J.S. 1990. The equilibrium partition function and base pair binding probabilities for RNA secondary structure. Biopolymers 29(6-7):1105-1119. https://doi.org/10.1002/bip.360290621 (PMID 1695107)
2. Lorenz, R. et al. 2011. ViennaRNA Package 2.0. Algorithms for Molecular Biology 6:26. https://doi.org/10.1186/1748-7188-6-26
3. Agarwal, V., Bell, G.W., Nam, J.W., Bartel, D.P. 2015. Predicting effective microRNA target sites in mammalian mRNAs. eLife 4:e05005. https://doi.org/10.7554/eLife.05005
4. Bernhart, S.H., Hofacker, I.L., Stadler, P.F. 2006. Local RNA base pairing probabilities in large sequences. Bioinformatics 22(5):614-615. https://doi.org/10.1093/bioinformatics/btk014 ; RNAplfold man page https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html
5. Turner, D.H., Mathews, D.H. NNDB: the nearest-neighbor parameter database (Turner 2004). https://rna.urmc.rochester.edu/NNDB/turner04/
