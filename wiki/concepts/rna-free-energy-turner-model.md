---
type: concept
title: "RNA free-energy model (Turner 2004 nearest-neighbor terms)"
tags: [rna, algorithm]
sources:
  - docs/Evidence/RNA-ENERGY-001-Evidence.md
  - docs/Evidence/RNA-HAIRPIN-001-Evidence.md
source_commit: 8346ce2d97d95f5b806caf203fd3d1dc19271cf5
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-energy-001-evidence
      evidence: "Test Unit: RNA-ENERGY-001 ... Title: Free Energy Calculation ... Area: RnaStructure"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-base-pairing
      source: rna-energy-001-evidence
      evidence: "Stacking energy parameters are assigned per adjacent base-pair stack; Watson-Crick stacking (16 entries) and G-U wobble stacking (20 entries) tables key off the WC {A-U,G-C} + G-U wobble pairing set"
      confidence: high
      status: current
---

# RNA free-energy model (Turner 2004 nearest-neighbor terms)

The **thermodynamic (energy) layer** of the RNA secondary-structure family (test unit
**RNA-ENERGY-001**): the individual **Turner 2004 nearest-neighbor free-energy terms** — stacking,
hairpin-loop initiation, terminal mismatch, dangling ends, special-hairpin totals, and AU/GU-end
penalties — that are summed to give the folding free energy `ΔG°` (kcal/mol at 37 °C). The records are
[[rna-energy-001-evidence]] (the aggregate energy layer) and [[rna-hairpin-001-evidence]] (a focused
deep-dive on the hairpin-loop + stem terms `CalculateHairpinLoopEnergy` / `CalculateStemEnergy`, with
NNDB worked Examples 1 and 2); [[test-unit-registry]] tracks the units, and
[[algorithm-validation-evidence]] describes the artifact pattern.

This is the **energy** layer, distinct from its two neighbours: the **chemistry** of *which* bases
pair ([[rna-base-pairing]] — Watson-Crick {A-U, G-C} + the single G-U wobble) and the **notation** a
fold is written in ([[rna-dot-bracket-notation]]). Stacking energies are assigned *to* the base pairs
of [[rna-base-pairing]], and these terms are what the folding units consume: the
[[pre-mirna-hairpin-detection|pre-miRNA hairpin]] `FreeEnergy` and the RNA-STRUCT-001 MFE folder both
sum them over a dot-bracket structure. Exposed methods on `RnaSecondaryStructure`:
`CalculateStackingEnergy` (one stack), `CalculateStemEnergy` (a stem's stacks), `CalculateHairpinLoopEnergy`
(loop destabilization), `CalculateMinimumFreeEnergy` (total MFE).

## 1. The nearest-neighbor model

Total folding free energy is a sum of local contributions:

```
ΔG°total = ΔG°initiation + Σ ΔG°stacking + Σ ΔG°loops
```

- **Stacking** — stabilizing (negative) energy from each pair of *adjacent* base pairs. **A single
  base pair has no stacking contribution (ΔG° = 0); stacking requires ≥2 adjacent pairs.**
- **Loops** — destabilizing (positive) energy from unpaired regions (hairpin, internal, bulge).
- **Initiation / ends** — cost of forming the first pair plus per-end AU/GU penalties.

All parameters are the **NNDB Turner 2004** tables at standard conditions **37 °C (310.15 K)** — a
defined condition, not an assumption. Precision is 2 decimal places (NNDB precision). Unknown /
non-canonical stacks contribute **0.0** (outside the NN model scope).

## 2. Stacking energies

**Watson-Crick stacks are all negative (stabilizing), and GC-rich stacks are the most stable:**

| Extreme | Stack (5'→3'/3'→5') | ΔG°37 (kcal/mol) |
|---------|---------------------|------------------|
| most stable | `GC/CG` | −3.42 |
| | `GG/CC` | −3.26 |
| least stable | `UA/AU` | −1.33 |
| | `AA/UU` | −0.93 |

**G-U wobble stacks are variable** — mostly negative (e.g. `GU/CG` −2.51, `CU/GG` −2.11) but **two are
destabilizing (positive)**: `UG/GU` **+0.30** and `GU/UG` **+1.29**. Two NNDB exceptions:

- **Note a:** `GG/UU` is set to **−0.5** for prediction.
- **Note b — special 3-stack context:** when `GU/UG` is flanked by GC and CG (`5'GGUC/3'CUGG`), the
  total for all three stacking interactions is **−4.12 kcal/mol**, *replacing* the individual sum
  (−1.53 + 1.29 + −1.53 = −1.77). Detected from the sequence.

## 3. Loop and end terms

- **Hairpin-loop initiation is positive and NOT monotonic in loop size** — an important invariant
  correction: 3→+5.4, 4→+5.6, 5→+5.7, **6→+5.4** (more stable than size 4), 7→+6.0, **8→+5.5**, 9→+6.4.
  For loop size **n > 9**: `ΔG°37(n) = ΔG°37(9) + 1.75·R·T·ln(n/9)` (Jacobson-Stockmayer extrapolation;
  R = 1.987 cal/mol·K, T = 310.15 K).
- **Special hairpins** (tri/tetra/hexaloops — e.g. **UNCG** and **GNRA** tetraloops) are poorly fit by
  the model, so NNDB supplies **total experimental energies keyed by closing pair + loop** (e.g.
  `CCUCGG` 2.5, the most stable UNCG) that **replace** the model calculation. GNRA tetraloops are more
  stable than random tetraloops.
- **All-C (poly-C) loops** carry an extra penalty from electrostatic repulsion: 3-nt flat **+1.5**;
  >3-nt **`0.3n + 1.6`**. Poly-C loops are less stable than equivalent poly-A loops.
- **Terminal mismatch** (96 NNDB entries) stabilizes the first unpaired pair adjacent to a helix end.
- **AU/GU-end penalty:** **+0.45 kcal/mol per AU end and per GU end** (same value for both), applied
  terminally and to inner AU/GU closures per the NNDB examples.

## 4. Worked oracles

- **GC stem (3 bp)** `GCG…CGC`: `GC/CG`(−3.42) + `CG/GC`(−2.36) = **−5.78 kcal/mol**.
- **NNDB hairpin example 1** `GGGAUAAAUCCC` (4-bp stem GC,GC,GC,AU + UAAA loop): stacking
  `GG/CC`(−3.26)+`GG/CC`(−3.26)+`GA/CU`(−2.35) = −8.87; hairpin init(4) 5.6 + terminal mismatch
  (AUAU −0.6) = 5.0; inner AU-end +0.45 → total **−3.42 kcal/mol**.
- **GGUC/CUGG special 3-stack:** **−4.12** (vs −1.77 without the special context).
- **NNDB hairpin example 1** (closing A-U, 6-nt loop A…A): loop +4.6 (init(6) 5.4 + terminal mismatch
  A·A −0.8) + helix −6.01 (3 stacks + AU-end) = total **−1.4** (RNA-HAIRPIN-001).
- **NNDB hairpin example 2** (closing A-U, 5-nt loop G…G): loop +4.1 (init(5) 5.7 − 0.8 tm − 0.8 GG
  first-mismatch bonus) + helix −6.01 = total **−1.9**. **3-nt loops get no first-mismatch term**;
  loops < 3 nt are prohibited; a stem of P pairs has P−1 stacks (P ≤ 1 ⇒ 0) (RNA-HAIRPIN-001).

## Invariants and edge cases

- **INV:** all Watson-Crick stacking energies are negative; hairpin-loop initiation is always positive;
  loop initiation is non-monotonic in size; GNRA tetraloops beat random tetraloops; poly-C < poly-A.
- **Empty sequence, single base pair, and poly-A (no structure) all give ΔG° = 0** — stacking needs ≥2
  adjacent pairs, so one pair (or none) contributes nothing.
- Null / empty input → 0 (graceful); invalid bases are an implementation choice (undefined); very long
  loops (>30) use the Jacobson-Stockmayer extrapolation.

## Scope and limitations

A [[scientific-rigor|research-grade]] reference for the **energy terms**, not a folder — assigning the
Turner ΔG° to a *given* structure is what these methods do; searching structure space for the minimum
is the RNA-STRUCT-001 MFE folder's job (which reuses these very terms). **No source contradictions** —
every parameter set is an exact match against **NNDB Turner 2004** (WC 16, GU 20 incl. notes a/b,
hairpin initiation 28, terminal mismatch 96, special tri/tetra/hexaloops 22, bonuses/penalties 6), and
Turner 2004, Mathews et al. 2004, Xia et al. 1998, and SantaLucia 1998 are mutually consistent. The
only recorded items are three **defined conditions** (37 °C standard state; 2-dp precision; unknown
stacks → 0.0), not assumptions.
