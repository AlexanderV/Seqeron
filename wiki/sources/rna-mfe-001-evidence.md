---
type: source
title: "Evidence: RNA-MFE-001 (Minimum Free Energy RNA folding — Zuker–Stiegler DP)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-MFE-001-Evidence.md
sources:
  - docs/Evidence/RNA-MFE-001-Evidence.md
source_commit: 18048c17104987bc34f26dbdfdc4c48f30fef2d2
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-MFE-001

The validation-evidence artifact for test unit **RNA-MFE-001** — **Minimum Free Energy (MFE) RNA
secondary-structure prediction**: the **Zuker–Stiegler dynamic-programming folder** that searches
structure space for the fold of lowest folding free energy (`CalculateMinimumFreeEnergy` /
`PredictStructure` on `RnaSecondaryStructure`). It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept is
[[rna-minimum-free-energy-folding]]. [[test-unit-registry]] tracks the unit.

This is the **folding / search** layer of the RNA secondary-structure family — the unit that
**consumes** the [[rna-free-energy-turner-model|Turner 2004 nearest-neighbor terms]] (summing them
over candidate structures) and **produces** a [[rna-dot-bracket-notation|dot-bracket]] structure as
output. Prior RNA ingests repeatedly flagged it as "the not-yet-ingested MFE folder that consumes the
Turner energy terms" (referred to generically as *RNA-STRUCT-001*); this Evidence file records the
unit under its own id **RNA-MFE-001**.

## What this file records

- **Algorithm:** MFE folding by the **Zuker–Stiegler (1981)** dynamic-programming decomposition with
  **Turner 2004** nearest-neighbor parameters. RNA structures are decomposed into distinct loop types —
  **hairpin loops**, **stacking regions** (stacked pairs), **bulge / interior loops**, and
  **multibranched (junction) loops** — and the total free energy is the sum of these loop
  contributions; the DP finds the decomposition of minimum total ΔG°.

- **Authoritative sources:**
  - **Zuker & Stiegler (1981)** *Nucleic Acids Research* 9(1):133–148 (rank 1, primary) — the original
    MFE DP; demonstrated on a 459-nt immunoglobulin mRNA fragment and *E. coli* 16S rRNA (polynomial,
    not exponential).
  - **Lorenz et al. (2011)** ViennaRNA Package 2.0, *Algorithms for Molecular Biology* 6:26 (rank 3) —
    ViennaRNA's MFE fold is "derived from the decomposition scheme … by Zuker and Stiegler [1981]";
    standard thermodynamic MFE folding runs in **O(n³)**.
  - **Ward, Datta, Wise, Mathews (2017)** *Nucleic Acids Research* 45(14):8541–8552 (rank 1) — the
    multiloop recurrences and complexity: the **standard affine (linear) multiloop Turner model runs in
    O(n³) time and O(n²) space**; the logarithmic multiloop model raises it to O(n⁴)/O(n³). Their title
    finding: *the simplest (affine) model is best*.
  - **NNDB Turner 2004** hairpin worked Examples 1 & 2, Watson-Crick helix example, and the
    stacking / loop-initiation / end-penalty parameter tables (rank 2) — the constants consumed by the
    folder, independently verified for the sibling unit RNA-HAIRPIN-001.
  - **Mathews et al. (2004)** PNAS 101:7287–7292 — primary parameter source.

- **DP matrices (Ward et al. 2017):** for every subsequence `(i,j)` compute **C(i,j)** = min free
  energy of the substructure *closed* by base pair `(i,j)` = `min(hairpin, interior/bulge over an inner
  pair, multiloop)`; multiloop fragment matrices **M** (≥1 component) and **M1** (exactly one
  component); and an exterior-loop matrix **F**. Seqeron implements the **standard affine multiloop
  model → O(n³) time / O(n²) space**.

- **Oracles (full single-hairpin MFE, from NNDB worked examples):**
  - `CalculateMinimumFreeEnergy("CACAAAAAAAUGUG")` = **−1.41 kcal/mol** (NNDB Example 1; NNDB rounds to
    −1.4). Unique optimal fold = 4-bp stem + 6-nt loop; `PredictStructure` yields dot-bracket
    `((((......))))` with pairs C-G / A-U / C-G / A-U.
  - `CalculateMinimumFreeEnergy("CACAGAAAGUGUG")` = **−1.91 kcal/mol** (NNDB Example 2, GG first
    mismatch; NNDB rounds to −1.9).

- **Boundary / corner cases:**
  - **Empty / null** sequence → **MFE = 0** (no structure).
  - **Unfoldable** input (homopolymer `AAAAAAAA`, no complementary bases) → **0** — the empty open-chain
    structure with ΔG° = 0 is always available, so the optimum is never positive.
  - Sequence **shorter than `minLoopSize + 2`** (e.g. `GCGC`, length < 5) → **0** — a hairpin must
    enclose ≥ **3** unpaired nt; pairs `(i,j)` with `j − i − 1 < 3` cannot close a hairpin.
  - **Intramolecular ⇒ no helix-initiation constant** — the bimolecular helix-initiation term is *not*
    added for a unimolecular fold (NNDB hairpin-example-1 note).

- **Invariants (property-based):**
  - **INV-01:** MFE ≤ 0 for every input (the 0-energy empty structure is always in the search set).
  - **INV-02:** monotonic non-increase under suffix extension — `MFE(s) ≤ MFE(prefix of s)`; extending
    a sequence only adds folding options.
  - **INV-03:** the optimized DP and the classic O(n³) baseline return identical scores under the
    comparable simplified pair-energy model (benchmark equality assertion across all lengths).

## Deviations and assumptions

Two documented items (neither invents a source-defined parameter):

1. **Multiloop per-unpaired cost `c = 0`** — the affine multiloop model `a + c·helices` uses offset
   `a = 9.25` and helix term `c = −0.63` (NNDB Turner 2004 multibranch parameters) with the
   per-unpaired-nucleotide coefficient set to **0** (`ML_unpaired = 0`), a documented simplification of
   the same affine family (Ward et al. 2017 confirm the affine model is the standard O(n³) choice). Does
   not affect the cited single-hairpin oracles (no multiloop present) but can shift the optimum for
   sequences that fold into multiloops. Marked a simplification, not an invented constant.
2. **Rounding** — NNDB tabulates final ΔG°37 to one decimal; the implementation rounds to two
   (`Math.Round(.., 2)`). Tests assert the exact two-decimal arithmetic sum with `.Within(1e-9)` **and**
   that `Math.Round(mfe, 1)` equals the one-decimal NNDB total, so no source-defined parameter changes.

**No source contradictions** — Zuker & Stiegler 1981, Lorenz 2011 (ViennaRNA), Ward 2017, and NNDB
Turner 2004 / Mathews 2004 are mutually consistent on the DP decomposition, the O(n³)/O(n²) affine
complexity, and the worked-example energies.
