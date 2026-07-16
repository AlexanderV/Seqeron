---
type: concept
title: "RNA minimum-free-energy folding (Zuker–Stiegler dynamic programming)"
tags: [rna, algorithm]
mcp_tools:
  - minimum_free_energy
sources:
  - docs/Evidence/RNA-MFE-001-Evidence.md
  - docs/Evidence/RNA-PARTITION-001-Evidence.md
  - docs/Evidence/RNA-STRUCT-001-Evidence.md
  - docs/algorithms/RnaStructure/Minimum_Free_Energy.md
source_commit: 2a42d224e78fd25e3f3a57b67a817e628ba2e33c
created: 2026-07-10
updated: 2026-07-16
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-mfe-001-evidence
      evidence: "Test Unit ID: RNA-MFE-001 ... Algorithm: Minimum Free Energy (MFE) RNA secondary structure prediction (Zuker–Stiegler dynamic programming with Turner 2004 nearest-neighbor parameters)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-free-energy-turner-model
      source: rna-mfe-001-evidence
      evidence: "ViennaRNA's MFE folding is 'derived from the decomposition scheme as described by Zuker and Stiegler'; the total free energy is the sum of Turner 2004 loop/stacking contributions the folder minimises over structure space (the folder CONSUMES the Turner NN terms)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-dot-bracket-notation
      source: rna-mfe-001-evidence
      evidence: "PredictStructure('CACAAAAAAAUGUG') yields dot-bracket ((((......)))) — the MFE folder produces a dot-bracket structure as output"
      confidence: high
      status: current
---

# RNA minimum-free-energy folding (Zuker–Stiegler dynamic programming)

The **folding / search layer** of the RNA secondary-structure family (test unit **RNA-MFE-001**):
given an RNA sequence, search structure space for the fold of **minimum total free energy** and return
both its ΔG° and its structure. This is the unit that **consumes** the
[[rna-free-energy-turner-model|Turner 2004 nearest-neighbor free-energy terms]] — summing them over
candidate structures — and **produces** a [[rna-dot-bracket-notation|dot-bracket]] structure as output.
The record is [[rna-mfe-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

Prior RNA ingests, before the RNA-STRUCT-001 artifact was ingested, flagged this as "the not-yet-ingested
MFE folder that consumes the Turner terms" and referred to it generically as *RNA-STRUCT-001* (the id the
[[pre-mirna-hairpin-detection]] opt-in `AssessHairpinByMfe` path and several source pages cite). That
alias assumption is now **superseded**: RNA-STRUCT-001 is a **distinct sibling unit** —
[[rna-secondary-structure-prediction|the top-level Nussinov base-pair-maximization predictor]]
(`Predict` / `PredictWithConstraints` / `ToDotBracket` / `FromDotBracket`) — while this physical
kcal/mol MFE folder is recorded under its own id **RNA-MFE-001**. The two are separate test units in the
`RnaStructure` area that **share** the same `RnaSecondaryStructure` Zuker MFE machinery (RNA-STRUCT-001's
deviation D5 added `CalculateMfeStructure`/`PredictStructureMfe`, the traceback partner of this unit's
`PredictStructure`), not two names for one unit.

## 1. What the algorithm does

MFE folding is solved by **dynamic programming** — the **Zuker & Stiegler (1981)** decomposition, the
same scheme ViennaRNA's MFE fold is derived from. Exposed on `RnaSecondaryStructure`:

- `CalculateMinimumFreeEnergy(seq)` — the **exact DP optimum** ΔG°37 (kcal/mol) of the implemented
  energy model.
- `PredictStructure(seq)` — a concrete fold as a dot-bracket structure + its base pairs, produced by a
  **greedy stem-loop assembly** (NOT the DP traceback), so its reported energy can be **above** the true
  MFE; for the optimal score use `CalculateMinimumFreeEnergy` (see §6).

The structure is **decomposed into distinct loop types**, and the total free energy is the *sum* of
their local contributions:

- **hairpin loops**, **stacking regions** (adjacent stacked pairs), **bulge / interior loops**, and
  **multibranched (junction) loops**.

The DP finds the decomposition of minimum total ΔG° — the Turner terms of
[[rna-free-energy-turner-model]] supply each contribution.

## 2. DP matrices and complexity (Ward et al. 2017)

For every subsequence `(i, j)` the folder computes:

- **C(i,j)** — min free energy of the substructure *closed* by base pair `(i, j)`, decomposed as
  `C(i,j) = min(hairpin, interior/bulge over an inner pair, multiloop)` ("a structure enclosed in a
  base pair is either a hairpin loop, delimited by an interior loop, or branches in a multiloop").
- **M / M1** — multiloop fragment matrices: `M` = a fragment with **one or more** components, `M1` = a
  fragment with **exactly one** component.
- **F** — the exterior-loop matrix.

Seqeron implements the **standard affine (linear) multiloop model** → **O(n³) time, O(n²) space**. The
alternative logarithmic multiloop model would raise this to O(n⁴)/O(n³) *and* change the optimum; Ward
et al. 2017 find the simplest (affine) model is best. The algorithm is polynomial, not exponential
(Zuker & Stiegler folded a 459-nt mRNA fragment and 16S rRNA fragments).

## 3. Worked oracles

Both from NNDB Turner 2004 worked examples, each constructed so its unique optimal fold is a single
hairpin:

| Sequence | Optimal structure | MFE ΔG°37 |
|----------|-------------------|-----------|
| `CACAAAAAAAUGUG` (len 14) | `((((......))))` — 4-bp stem C-G/A-U/C-G/A-U + 6-nt loop | **−1.41** (NNDB rounds −1.4) |
| `CACAGAAAGUGUG` (len 13) | 4-bp stem + 5-nt loop, GG first mismatch | **−1.91** (NNDB rounds −1.9) |

Example 1 term breakdown: stacks −2.11, −2.24, −2.11; AU end +0.45; terminal mismatch (AU·AA) −0.8;
hairpin initiation(6) +5.4. (Helix −6.01 + loop +4.6 = −1.41.)

## 4. Invariants and edge cases

- **INV-01 (MFE ≤ 0):** the empty open-chain structure with ΔG° = 0 is always in the search set, so the
  optimum is never positive.
- **INV-02 (suffix monotonicity):** `MFE(s) ≤ MFE(prefix of s)` — extending a sequence only adds
  folding options (monotonic non-increase under suffix extension).
- **INV-03 (engine agreement):** the optimized DP and the classic O(n³) baseline return identical
  scores under the comparable simplified pair-energy model (benchmark equality across all lengths).
- **No structure ⇒ 0:** empty / null input, an unfoldable homopolymer (`AAAAAAAA` — no complementary
  bases), and any sequence **shorter than `minLoopSize + 2`** (e.g. `GCGC`, length < 5) all return
  **MFE = 0**. A hairpin must enclose **≥ 3** unpaired nt; pairs `(i, j)` with `j − i − 1 < 3` cannot
  close a hairpin.
- **Intramolecular ⇒ no helix-initiation constant** — the bimolecular helix-initiation term is *not*
  added for a unimolecular fold (NNDB hairpin-example-1 note).

## 5. Scope, assumptions, and relationships

A [[scientific-rigor|research-grade]] MFE folder — it **searches** structure space for the *single*
optimal fold, whereas the [[rna-partition-function-mccaskill|McCaskill partition function]] is its
Boltzmann-weighted **ensemble** counterpart (a probability distribution over all structures + per-pair
probabilities rather than one optimum). By contrast [[rna-free-energy-turner-model]] only **scores** a
*given* structure with the Turner terms this folder consumes, and [[rna-dot-bracket-notation]] only
**parses/validates** the notation this folder emits. Both this folder and the McCaskill ensemble
search only **nested (pseudoknot-free)** structure space; the crossing-helix extension —
**pseudoknots** — is [[rna-pseudoknot-prediction]] (RNA-PKPREDICT-001), which **depends on** this
folder both as its fallback baseline (`FreeEnergy ≤ CalculateMfeStructure().FreeEnergy`) and to fold
its internal loops.
It builds on the [[rna-base-pairing]] Watson-Crick + G-U wobble pairing rule. The
[[pre-mirna-hairpin-detection|pre-miRNA hairpin detector]]'s opt-in `AssessHairpinByMfe` /
`FindPreMiRnaHairpinsByMfe` paths fold with this engine and read the hairpin from the real MFE
structure.

Two documented items (neither invents a source parameter):

1. **Multiloop per-unpaired cost `c = 0`** — the affine multiloop model `a + c·helices` uses
   `a = 9.25`, helix term `c = −0.63` (NNDB Turner 2004 multibranch) with the per-unpaired-base
   coefficient set to 0 (`ML_unpaired = 0`); a documented simplification of the same affine family. Does
   not affect the single-hairpin oracles (no multiloop) but can shift the optimum for multiloop folds.
2. **Rounding** — implementation rounds to two decimals; tests assert the exact two-decimal arithmetic
   sum `.Within(1e-9)` and that `Math.Round(mfe, 1)` equals the one-decimal NNDB total.

**No source contradictions** — Zuker & Stiegler 1981, Lorenz 2011 (ViennaRNA), Ward 2017, and NNDB
Turner 2004 / Mathews 2004 agree on the DP decomposition, the O(n³)/O(n²) affine complexity, and the
worked-example energies.

## 6. Implementation surface (RNA-MFE-001 primary spec)

The canonical algorithm spec is `docs/algorithms/RnaStructure/Minimum_Free_Energy.md`; the code lives in
`RnaSecondaryStructure` in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`
(Implementation Status **Simplified**). Three entry points — two public, one internal comparator:

- **`CalculateMinimumFreeEnergy(string rnaSequence, int minLoopSize = 3)` → `double`** — the Zuker-style
  O(n³) DP score (ΔG°37, kcal/mol, **rounded to 2 decimals**, ≤ 0). `minLoopSize` is **clamped up to 3**
  if smaller. It fills `C(i,j)` / `WM(i,j)` by increasing span (hairpin, stacking, interior/bulge,
  multibranch), accumulates the exterior optimum `W(j)` with AU/GU end penalties and dangling ends, and
  returns `W(n−1)`. This is the **exact global optimum of the affine-multiloop / no-coaxial-stacking
  energy model** — not a heuristic.
- **`PredictStructure(string, int minLoopSize = 3, int minStemLength = 3, int maxLoopSize = 10)` →
  `SecondaryStructure`** — a **greedy** heuristic that selects non-overlapping stem-loops by energy;
  O(n²) per stem-loop scan, O(n) space. It may return a structure whose energy is **above** the true
  MFE. The returned record carries `Sequence`, `DotBracket`, `BasePairs`, `StemLoops`, `Pseudoknots`,
  and `MinimumFreeEnergy`.
- **`CalculateMinimumFreeEnergyClassic(string, int)`** (internal) — an O(n³) baseline using a
  **simplified per-pair (Nussinov-style) energy model**, retained only as the **timing comparator** in
  the benchmark. It does **not** use the Turner model and is **not expected to match** the Turner-model
  score (this is the "classic baseline" of INV-03; the doc-header flags the whole classic path as a
  *baseline/reference* comparator, not full thermodynamic fidelity).

Complexity / memory details the algorithm spec pins down:

- The interior/bulge search is bounded by **`MAXLOOP = 30`** (sliding lower bounds) → a constant window,
  so it is **O(1) per `(i,j)`**; the **multiloop split is the O(n³) term** (§4.3).
- The `n×n` DP matrices are a **flat `double[]` buffer indexed `i*n + j`, rented from
  `ArrayPool<double>`** → O(n²) space.
- The **GGUC/CUGG special 3-stack context** (NNDB note b) is detected from the sequence and scored as a
  single **−4.12** term (consistent with [[rna-free-energy-turner-model]] §2).
- **Suffix-tree reuse was evaluated and is NOT used** — MFE is a thermodynamic-scoring DP, not an
  exact-substring search, so the suffix tree offers no benefit for energy minimization (recorded per the
  "reuse existing infrastructure" rule).

**Preconditions / edge behaviour (no exceptions thrown):** `null`/empty → 0; a sequence shorter than
`minLoopSize + 2` → 0 (cannot enclose a 3-nt hairpin); input is upper-cased and 0-based; only A/C/G/U
pair (WC + G-U wobble) and unrecognized characters simply never pair.

**Performance baseline** (full Turner model, `RnaSecondaryStructure_MFE_Benchmark`, median of 3, Release,
Apple-silicon, 2026-06-14) — confirms the §2 cubic growth:

| Length n | 100 | 200 | 300 | 500 | 1000 |
|----------|-----|-----|-----|-----|------|
| Time (ms) | ~25 | ~83 | ~214 | ~639 | ~3155 |

**Deviations** the spec records (neither invents a source parameter): (1) multiloop per-unpaired cost
`= 0` (ASM-04) — shifts only multiloop optima, single-hairpin/stack/interior optima unaffected; (2)
**2-decimal rounding vs NNDB 1-decimal** — display precision only (tests assert `.Within(1e-2)` vs the
NNDB total and `.Within(1e-9)` vs the exact two-decimal sum).
