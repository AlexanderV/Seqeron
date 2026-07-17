---
type: concept
title: "Turner-2004 McCaskill partition function (ensemble accessibility & TargetScan SA)"
tags: [rna, algorithm]
sources:
  - docs/algorithms/RnaStructure/Turner_McCaskill_Partition_Function.md
source_commit: 18007bb74684cbe34dac3ee21116ec1f15bf37fb
created: 2026-07-17
updated: 2026-07-17
---

# Turner-2004 McCaskill partition function (ensemble accessibility & TargetScan SA)

The **full-thermodynamic** McCaskill (1990) partition function evaluated on the **same Turner-2004
nearest-neighbour energy model** the library's Zuker–Stiegler [[rna-minimum-free-energy-folding|MFE
folder]] uses — not the base-pair-counting teaching model. Instead of taking the minimum free energy
over all pseudoknot-free structures, it Boltzmann-sums `exp(−E(S)/RT)` over the whole ensemble to
yield the partition function `Z`, the equilibrium base-pair probabilities `P(i,j)`, the **per-base
unpaired probabilities** `p_unpaired(i) = 1 − Σ_j P(i,j)`, the **ensemble free energy** `−RT·ln Z`,
and the **joint region-unpaired probability** (RNAplfold-style accessibility). That accessibility
feeds the TargetScan context++ **SA** (structural accessibility) feature used by the miRNA
target-site scorer. Test unit **MIRNA-TARGET-001** (SA wiring) / a new RNA-STRUCT capability,
implemented in `RnaSecondaryStructure.cs` in `Seqeron.Genomics.Analysis`. See
[[test-unit-registry]] for how the unit is tracked and [[algorithm-validation-evidence]] for the
artifact pattern.

## Why this is a distinct unit (not the simplified McCaskill page)

This is the **full-Turner sibling** of the simplified base-pair-counting
[[rna-partition-function-mccaskill]] (RNA-PARTITION-001, `CalculatePartitionFunction`). That page
documents the same McCaskill recurrence over a **fixed per-pair energy `E_bp`** (the Freiburg
teaching model) and explicitly lists a *full Turner-parameter partition function* as **not
implemented**. This unit is that missing engine — same pseudoknot-free structure space and Boltzmann
formalism, but a genuinely different implementation:

| Aspect | Turner-2004 McCaskill (**this**) | Simplified McCaskill ([[rna-partition-function-mccaskill]]) |
|---|---|---|
| Energy model | full **Turner-2004** NN (stacking, hairpin/internal/bulge init + mismatches, dangles, terminal-AU, multibranch `a/b/c`) | fixed `E_bp` per pair (teaching model) |
| Entry points | `CalculateUnpairedProbabilities`, `CalculateRegionUnpairedProbability` | `CalculatePartitionFunction`, `CalculateStructureProbability` |
| Outputs | `Z`, `P(i,j)`, **`p_unpaired`**, **`ΔG_ensemble`**, **region accessibility** | `Z`, `P(i,j)` |
| Pair-prob method | **constrained re-folds** (`Z_require(i,j)/Z`), no outside recursion | full **outside** recursion `Q^b·O/Z` |
| Complexity | `O(n⁵)` worst-case bpp, `O(n³)` for `Z`/`p_unpaired`/accessibility | `O(n³)` |
| Application | miRNA **SA** feature (TargetScan context++) | structure counting / ensemble teaching oracles |

Because it exposes distinct outputs, distinct entry points, a distinct implementation strategy, and
a distinct application (SA), it is its own page rather than an implementation surface on the
simplified concept.

## Boltzmann-weighted nearest-neighbour recurrences

The inside recursion mirrors the MFE Zuker recurrence under the substitution `min → Σ` and
`+ΔG → × exp(−ΔG/RT)`, reusing the **same** Turner-2004 loop-energy functions as the MFE `FillDp`:

- `Vexp(i,j)` — the sub-ensemble in which `i·j` are paired: Boltzmann-sums hairpin / stacking /
  internal / bulge / multibranch closings.
- `WMexp` — the multibranch (multiloop) region.
- `Wexp` — the external loop; the total partition function is `Z = Wexp(0, n−1)`.

Constants: `β = 1/(RT)`, `RT = R·T` with `R ≈ 1.987×10⁻³ kcal/(mol·K)`, default `T = 310.15 K`
(37 °C, the Turner/NNDB standard) so `RT = 0.61626805 kcal/mol`. Definitions:

```text
Z            = Σ_{s∈Ω} exp(−β·E(s))          (Ω = pseudoknot-free structures)
p(s)         = exp(−β·E(s)) / Z
P(i,j)       = (Σ_{s ∋ (i,j)} exp(−β·E(s))) / Z
p_unpaired(i)= 1 − Σ_j P(i,j)
ΔG_ensemble  = −RT·ln Z
P([a..b] unpaired) = Z_open(a,b) / Z
```

## Constrained-refold implementation (no outside recursion)

Rather than a bespoke outside recursion (the route the simplified page uses for `P(i,j)`), this
unit computes both marginals by **constrained re-folds of the same inside DP**:

- `p_unpaired(i) = Z_forbid(i) / Z` — re-fill the inside matrices while forbidding position `i` from
  pairing.
- `P(i,j) = Z_require(i,j) / Z` — re-fill requiring the pair `(i,j)` with proper non-crossing
  nesting.

This guarantees the normalisation `Σ_j P(i,j) + p_unpaired(i) = 1` to floating-point precision
(verified to ~1e-15) exactly, at the cost of one `O(n³)` fill per base (for `p_unpaired`) and per
candidate pair (for `P(i,j)`) — hence the `O(n⁵)` worst case for the full base-pair-probability
matrix. `Z`, `p_unpaired`, and region accessibility remain `O(n³)`; space is `O(n²)`.

## Region accessibility and the miRNA SA feature

`CalculateRegionUnpairedProbability(seq, windowEnd, windowLength, …)` returns the **joint**
probability that a whole window `[a..b]` is unpaired, `Z_open(a,b)/Z` — RNAplfold's `_lunp` row `i`,
column `L` reports `P([i−L+1 .. i] unpaired)`. `MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus` wires
this into the TargetScan context++ **SA** term: for a seed-matched site it reads the accessibility
of the **14-nt window ending 7 nt 3′ of the seed-match start** (7mer-A1 decrements the start by 1),
`log10`-transforms it, and min-max scales by the verbatim Agarwal-2015 SA coefficients. The local
fold uses a single `W = 80` context window (bases pair only within ±`L = 40`) — RNAplfold reference
parameters `-L 40 -W 80 -u 20`. When the SA window does not fit the UTR, SA contributes 0 and is
reported in `OmittedFeatures` (matching the perl reference's missing-`_lunp` → 0). This is the SA
input consumed by [[mirna-target-site-prediction]]'s context++ scorer.

## Contract

`CalculateUnpairedProbabilities(rnaSequence, minLoopSize = 3, temperature = 310.15)` returns an
`UnpairedProbabilityResult`:

| Field | Type | Meaning |
|---|---|---|
| `PartitionFunction` | `double` | `Z` |
| `BasePairProbabilities` | `IReadOnlyDictionary<(int,int),double>` | `P(i,j)` for 0-based `i<j` (only `P>0`) |
| `UnpairedProbabilities` | `IReadOnlyList<double>` | `p_unpaired(i)`, one per position |
| `EnsembleFreeEnergy` | `double` | `−RT·ln Z` (kcal/mol) |

`CalculateRegionUnpairedProbability(...)` returns a single `double` — `P(window unpaired) ∈ [0,1]`.
Coordinates are **0-based inclusive**. Inputs: `minLoopSize` clamped to ≥ 3 (NNDB steric floor);
`temperature` in absolute K, must be `> 0`; the accessibility window (0-based inclusive end +
length > 0) must fit the sequence. DNA is accepted (`T` read as `U`), upper-cased internally.

## Invariants (test oracles)

- **INV-01** `Z ≥ 1` (`> 0` always) — the open chain (`E = 0`) always contributes `exp(0) = 1`.
- **INV-02** `P(i,j) ∈ [0,1]` — a ratio of a non-negative sub-sum to `Z`.
- **INV-03** `Σ_j P(i,j) ≤ 1` — at most one pair involves position `i` in any nested structure.
- **INV-04** `p_unpaired(i) = 1 − Σ_j P(i,j) ∈ [0,1]`.
- **INV-05** `ΔG_ensemble = −RT·ln Z ≤ MFE` — the MFE structure is one term in `Z`.
- **INV-06** region accessibility `Z_open/Z ∈ [0,1]` (`Z_open` sums a subset of `Z`).

## Worked example (`GAAAC`)

`GAAAC` can adopt only the open chain or the single hairpin closed by G(0)–C(4) over the 3-nt loop
`AAA` (Turner hairpin energy 5.4 kcal/mol; G–C ⇒ no terminal-AU penalty), with `RT = 0.61626805`:

- `Z = 1 + exp(−5.4/RT) = 1.0001565052764922`
- `P(0,4) = exp(−5.4/RT)/Z = 0.00015648078642340854`
- `p_unpaired(0) = p_unpaired(4) = 0.9998435192135765`; `p_unpaired(2) = 1`
- `ΔG_ensemble ≈ −9.644×10⁻⁵ ≤ MFE = 0`
- `CalculateRegionUnpairedProbability("GAAAC", 4, 5) = 1/Z` (only the open chain leaves the whole
  window unpaired).

## Edge cases and scope

null / empty → `Z = 1`, no pairs, `ΔG = 0` (only the empty structure); too short to pair → `Z = 1`,
all bases unpaired; `temperature ≤ 0` → `ArgumentOutOfRangeException`; an accessibility window that
does not fit → `ArgumentOutOfRangeException` (a window in a too-short sequence returns `1.0`); the
SA window off the UTR end → SA = 0, listed in `OmittedFeatures`.

**Modelling assumptions / limitations** — a [[scientific-rigor|research-grade]] ensemble method.
Nearest-neighbour additivity of the Turner-2004 model (ensemble probabilities inherit any model
error); **pseudoknot-free** only (crossing pairs are excluded from `Z` and accessibility — no NN
ensemble model includes them, see [[rna-pseudoknot-detection]]); probabilities are
**temperature-specific** (a different `T` re-weights the ensemble). The per-pair bpp refold is
`O(n⁵)`, intended for the short windows SA uses, not genome-scale folding; the SA local fold is a
single `W = 80` window rather than RNAplfold's full sliding-window average, so values can differ
slightly near long-range context. PCT (conservation) is not part of this capability.

## Relationships

It is the **`alternative_to` / ensemble counterpart** of the single-optimum
[[rna-minimum-free-energy-folding]] MFE folder over the *same* Turner-2004 model, and the
**full-thermodynamic sibling** of the simplified [[rna-partition-function-mccaskill]]. It
Boltzmann-weights the [[rna-free-energy-turner-model]] nearest-neighbour free-energy terms, builds
on the [[rna-base-pairing]] Watson-Crick + G-U wobble rule, and shares the pseudoknot-free nesting
convention with the whole RNA secondary-structure family. Its region-accessibility output is the SA
input that [[mirna-target-site-prediction]]'s TargetScan context++ scorer consumes.
