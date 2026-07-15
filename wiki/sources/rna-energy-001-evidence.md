---
type: source
title: "Evidence: RNA-ENERGY-001 (RNA free-energy calculation — Turner 2004 nearest-neighbor terms)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-ENERGY-001-Evidence.md
sources:
  - docs/Evidence/RNA-ENERGY-001-Evidence.md
source_commit: 79389a8898bca00ca6eabf4fef4b51f3ea73fae6
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-ENERGY-001

The validation-evidence artifact for test unit **RNA-ENERGY-001** — **Free Energy Calculation**, the
**thermodynamic (energy) layer of the RNA secondary-structure family** (`CalculateStemEnergy`,
`CalculateHairpinLoopEnergy`, `CalculateStackingEnergy`, `CalculateMinimumFreeEnergy` on
`RnaSecondaryStructure`). It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept is
[[rna-free-energy-turner-model]]. [[test-unit-registry]] tracks the unit.

The energies are assigned to adjacent RNA base-pair stacks, so the unit sits directly on the
base-pairing chemistry primitive [[rna-base-pairing]] (Watson-Crick {A-U, G-C} + the single G-U
wobble), and it is the energy the folding units consume — the [[pre-mirna-hairpin-detection|pre-miRNA
hairpin]] `FreeEnergy` and the RNA-STRUCT-001 MFE folder both sum these Turner terms over a
[[rna-dot-bracket-notation|dot-bracket]] structure.

## What this file records

- **Nearest-neighbor (NN) model** — total folding free energy is
  `ΔG°total = ΔG°initiation + Σ ΔG°stacking + Σ ΔG°loops`: an initiation cost for the first pair,
  stabilizing (negative) stacking from each adjacent base-pair stack, and destabilizing (positive)
  loop terms for unpaired regions (hairpin / internal / bulge).

- **Authoritative sources (all reference-implementation / curated-database / primary literature):**
  - **NNDB — Nearest Neighbor Database (Turner 2004)**, rna.urmc.rochester.edu/NNDB/turner04/ — the
    definitive parameter tables at 37 °C: WC stacking (16 entries), G-U stacking (20, incl. notes a/b),
    hairpin-loop initiation (`loop.txt`), terminal mismatch (`tm-parameters.html`, 96 entries), special
    tri/tetra/hexaloop total energies (`tloop.txt`/`triloop.txt`/`hexaloop.txt`, 22 entries), hairpin
    bonuses/penalties (`hairpin-mismatch-parameters.html`).
  - **Turner (2004)** *Thermodynamics of RNA* (RNA World, 3rd Ed.); **Mathews et al. (2004)** PNAS
    101(19):7287-7292; **Xia et al. (1998)** Biochemistry 37(42):14719-35; **SantaLucia (1998)** PNAS
    95(4):1460-5; **Heus & Pardi (1991)** Science 253:191-194; Wikipedia nucleic-acid thermodynamics /
    secondary structure.
  - **Reference implementations:** ViennaRNA (Turner parameters), MFOLD, RNAstructure (Mathews Lab).

- **Parameter tables (exact values worth pinning):**
  - **Watson-Crick stacking ΔG°37 (kcal/mol), all negative:** most stable `GC/CG` −3.42 and `GG/CC`
    −3.26 (GC-rich stacks); least stable `AA/UU` −0.93, `AU/UA` −1.10, `UA/AU` −1.33.
  - **G-U wobble stacking is variable** — mostly negative (e.g. `GU/CG` −2.51) but **two are
    destabilizing (positive)**: `UG/GU` +0.30 and `GU/UG` +1.29. Note a: `GG/UU` set to −0.5 for
    prediction. **Note b special context:** `5'GGUC/3'CUGG` (GU/UG flanked by GC and CG) — the total
    for all 3 stacking interactions is **−4.12 kcal/mol**, replacing the individual sum (−1.77).
  - **Per AU end and per GU end:** +0.45 kcal/mol (terminal / inner AU-GU penalty; same value for both).
  - **Hairpin-loop initiation is positive and NOT monotonic** in loop size: 3→+5.4, 4→+5.6, 5→+5.7,
    6→+5.4, 7→+6.0, 8→+5.5, 9→+6.4 (size 6 is more stable than size 4). For n > 9,
    `ΔG°37(n) = ΔG°37(9) + 1.75·R·T·ln(n/9)` (Jacobson-Stockmayer; R = 1.987 cal/mol·K, T = 310.15 K).
  - **Special hairpins** (tri/tetra/hexaloops, e.g. UNCG/GNRA tetraloops) use **NNDB total experimental
    energies keyed by closing pair + loop** (e.g. `CCUCGG` 2.5, most stable UNCG) that **replace** the
    model calculation.
  - **All-C loop penalty** (electrostatic repulsion): 3-nt flat +1.5; >3-nt `0.3n + 1.6`.

- **Datasets / oracles:**
  - **GC stem (3 bp)** `GCG…CGC`: `GC/CG`(−3.42) + `CG/GC`(−2.36) = **−5.78**.
  - **NNDB hairpin example 1** `GGGAUAAAUCCC` (4-bp stem GC,GC,GC,AU + UAAA loop): stacking
    `GG/CC`(−3.26)+`GG/CC`(−3.26)+`GA/CU`(−2.35) = −8.87; hairpin init(4) 5.6 + terminal-mismatch
    (AUAU −0.6) = 5.0; inner AU-end +0.45 → total **−3.42 kcal/mol**.
  - **GGUC/CUGG special 3-stack:** total **−4.12** (vs −1.77 without the special context).

- **Edge cases / invariants:** empty sequence, single base pair, and poly-A (no structure) all give
  **ΔG° = 0** (stacking needs ≥2 adjacent pairs — a single pair contributes 0); all-C loop → penalty;
  very long loop (>30) → Jacobson-Stockmayer extrapolation. Invariants: WC stacking all negative;
  hairpin-loop initiation always positive; GNRA tetraloops more stable than random; poly-C less stable
  than poly-A.

## Deviations and assumptions

Recorded as **not-assumptions** — every parameter is sourced from NNDB Turner 2004, so the file states
"No assumptions remain":

- **Standard conditions 37 °C (310.15 K)** — defined by Turner 2004, not chosen.
- **Precision 2 decimal places** — follows NNDB parameter precision.
- **Unknown / non-canonical stacking pairs contribute 0.0 kcal/mol** — outside the NN model scope.

**No source contradictions** — all parameter sets verified as exact matches against NNDB Turner 2004
(WC 16, GU 20, hairpin initiation 28, terminal mismatch 96, special hairpins 22, bonuses/penalties 6).
