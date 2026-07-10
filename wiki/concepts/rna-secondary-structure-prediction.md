---
type: concept
title: "RNA secondary-structure prediction (Nussinov base-pair maximization + constraints)"
tags: [rna, algorithm]
sources:
  - docs/Evidence/RNA-STRUCT-001-Evidence.md
source_commit: bb82b7ec80bbbf5750e53616ccc60df7b45c010c
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-struct-001-evidence
      evidence: "Test Unit: RNA-STRUCT-001 ... Title: Secondary Structure Prediction ... Area: RnaStructure"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:rna-minimum-free-energy-folding
      source: rna-struct-001-evidence
      evidence: "Both fold a sequence to a structure but with different objectives: RNA-STRUCT-001 is a Nussinov base-pair-maximizing DP with weighted pair scores (WC −2.0, Wobble −1.0) that 'Maximizes weighted pair count (not thermodynamic MFE)' and whose 'Results indicate relative stability, not physical energy units', whereas RNA-MFE-001 is the physical Turner-2004 kcal/mol Zuker–Stiegler MFE folder; D5 added CalculateMfeStructure/PredictStructureMfe so both units share the same V/W/WM matrices"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-base-pairing
      source: rna-struct-001-evidence
      evidence: "Must Test 1 'Base pairing correctness (Watson-Crick, Wobble)'; Base Pairing Invariants: Watson-Crick pairs A-U/U-A/G-C/C-G, Wobble G-U/U-G — the Nussinov DP maximizes over exactly these pairs"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-dot-bracket-notation
      source: rna-struct-001-evidence
      evidence: "Method coverage ToDotBracket(structure) (Notation — Output format) and FromDotBracket(notation) (Parse — Input format parsing); Must Test 4 'Dot-bracket notation correctness'; Simple Hairpin expected dot-bracket ((((....))))"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-free-energy-turner-model
      source: rna-struct-001-evidence
      evidence: "Stem-Loop Energy Model uses Turner 2004 nearest-neighbor stacking parameters from NNDB (terminal AU/GU +0.45, hairpin-loop initiation 3–30, special tri/tetra/hexaloop overrides, UU/GA −0.9 & GG −0.8 mismatch bonuses, all-C loop penalty) — §3.3 Energy Parameters (Turner 2004 — NNDB)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-pseudoknot-detection
      source: rna-struct-001-evidence
      evidence: "§6.3 Pseudoknot Detection: crossing (0,6),(3,9) → pseudoknot detected; nested (0,10),(2,8),(4,6) → none — the same i<k<j<l crossing test RNA-PSEUDOKNOT-001 scans (prediction is out of scope: 'No pseudoknot prediction — O(n³) DP inherently cannot find pseudoknots; correct by design')"
      confidence: high
      status: current
---

# RNA secondary-structure prediction (Nussinov base-pair maximization + constraints)

The **top-level structure-prediction umbrella** of the RNA secondary-structure family (test unit
**RNA-STRUCT-001**, area `RnaStructure`): given an RNA sequence, predict its secondary structure by a
**Nussinov base-pair-maximizing** dynamic program, optionally under **forced-pair constraints**, and
convert between the structure and its [[rna-dot-bracket-notation|dot-bracket]] string. The canonical
surface on `RnaSecondaryStructure` is `Predict(sequence)`; the record is [[rna-struct-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the artifact
pattern.

## 1. What distinguishes this unit — and its relationship to RNA-MFE-001

This is **not** the physical minimum-free-energy folder. It is a distinct sibling that solves the
"predict a structure from a sequence" problem with a **different objective**:

| Unit | Objective | Score | Units |
|------|-----------|-------|-------|
| **RNA-STRUCT-001** (this) | maximize **(weighted) base-pair count** — Nussinov | WC −2.0 / wobble −1.0 | *relative stability, not physical energy* |
| [[rna-minimum-free-energy-folding\|RNA-MFE-001]] | minimize **total free energy** — Zuker–Stiegler | Turner 2004 NN | kcal/mol (physical ΔG°) |

Both live on the same `RnaSecondaryStructure` class and **share the Zuker MFE machinery**: this unit's
**deviation D5** (2026-06-23) added `CalculateMfeStructure` / `PredictStructureMfe` (a Zuker–Stiegler
DP traceback over the *same* V/W/WM matrices as the scalar MFE), while the greedy `Predict` /
`PredictStructure` is kept as a fast path. So RNA-STRUCT-001 and RNA-MFE-001 are **two test units in one
area exercising overlapping code**, not two names for one unit — correcting an earlier assumption (made
before this artifact was ingested) that treated *RNA-STRUCT-001* as a generic alias of the MFE folder.

## 2. The algorithms

- **Nussinov & Jacobson (1980)** — the O(n³) time / O(n²) space DP that **maximizes the number of base
  pairs** over the [[rna-base-pairing|Watson-Crick + G-U wobble]] pairing set. Cannot represent
  pseudoknots (the recurrences are pseudoknot-free by construction — a correct-by-design boundary, not
  a bug).
- **Weighted-pair MFE DP** — the same Nussinov recursion with **weighted** pair scores (WC **−2.0**,
  wobble **−1.0**), maximizing weighted pair count. Its output signals **relative stability**, not
  kcal/mol; contrast the physical Turner ΔG° of RNA-MFE-001.
- **Stem-loop energy model** — Turner 2004 / NNDB nearest-neighbor terms (stacking, terminal AU/GU
  +0.45 per helix end, hairpin-loop initiation for sizes 3–30 with Jacobson–Stockmayer beyond, special
  tri/tetra/hexaloop total-energy overrides, first-mismatch bonuses UU/GA −0.9 & GG −0.8, all-C loop
  penalty). This is the [[rna-free-energy-turner-model]] term set, shared across the family.
- **Zuker traceback** (D5) — reconstructs the optimal structure whose energy equals
  `CalculateMinimumFreeEnergy` for the same input (asserted across hairpin / multi-stem / multiloop).

## 3. Method surface

`Predict(sequence)` (complete prediction), `PredictWithConstraints(seq, constraints)` (**forced base
pairs** — the Mathews et al. 2004 constrained-DP class), `ToDotBracket(structure)` (structure →
notation), `FromDotBracket(notation)` (notation → structure). The last two are the round-trip bridge to
[[rna-dot-bracket-notation]]; structural motifs classified are stem / hairpin loop / internal loop /
bulge / multi-loop / pseudoknot (the last **detected**, never predicted).

## 4. Worked oracles

| Sequence | Expected | Note |
|----------|----------|------|
| `GGGGAAAACCCC` | `((((....))))`, 4-bp stem + 4-nt loop, MFE < 0 | simple hairpin |
| `GCGCGAAACGCGC` | GNRA (GAAA) with G-C closing → standard model + GA first-mismatch bonus −0.9 | no NNDB special entry for GNRA/C-G |
| `GCGGAUUUAGCUCAGUUGG…GAAUUCGCA` (72 nt) | multiple stem-loops, cloverleaf, valid dot-bracket | tRNA-like |
| `AAAAAAAAAAAA` | no pairs, MFE = 0, all-dots | A cannot pair with A |

Pseudoknot **detection** oracles: crossing (0,6)+(3,9) → detected; non-crossing (0,5)+(1,4) → none;
nested (0,10)+(2,8)+(4,6) → none — the `i<k<j<l` crossing test shared with [[rna-pseudoknot-detection]].

## 5. Invariants and edge cases

- **Dot-bracket balance:** openers = closers; each base in ≤ 1 pair; selected stem-loops non-overlapping.
- **MFE sign ≤ 0** (Nussinov weighted-pair score); WC stacking always negative; loop initiation always
  positive (destabilizing).
- **Base pairing:** WC {A-U, U-A, G-C, C-G} + wobble {G-U, U-G} pair; all seven other combinations do
  not (shared [[rna-base-pairing]] rule).
- **No structure ⇒ 0:** empty `""` / null → empty structure, MFE = 0; too-short `"GC"` → no stem-loop;
  poly-A / poly-U → no structure. **Minimum hairpin** = 3-bp stem + 3-nt loop (steric floor). Case
  insensitive; wobble-only stems fold iff wobble enabled; invalid characters handled or rejected.

## 6. Scope, limitations, relationships

A [[scientific-rigor|research-grade]] predictor. Documented scope boundaries (all accepted, not bugs):
**no pseudoknot prediction** (the O(n³) DP cannot represent crossings — detection only, via the
[[rna-pseudoknot-detection]] primitive; the family's energy-based crossing predictor is
[[rna-pseudoknot-prediction]]); **single sequence only** (no comparative/covariance folding); **minimum
loop size 3** (NNDB steric constraint). Two open items are the **int21 (2,304-entry)** and **int22
(36,864-entry)** internal-loop lookup tables — too large for inline static data, so a generic
initiation + asymmetry + mismatch model is used instead (deviations D3/D4, BLOCKED pending an external
data file).

Within the family it is the **alternative_to** the physical [[rna-minimum-free-energy-folding|MFE
folder]] (same problem, base-pair-max vs energy-min objective), builds on [[rna-base-pairing]], converts
to/from [[rna-dot-bracket-notation]], consumes the [[rna-free-energy-turner-model]] stem-loop terms, and
shares its crossing-detection test with [[rna-pseudoknot-detection]]. **No source contradictions** —
Nussinov & Jacobson 1980, Zuker & Stiegler 1981, MIT 6.047, Turner 2004 / NNDB, and Mathews 2004 are
mutually consistent; the only reconciliation is the RNA-STRUCT-001 ≠ RNA-MFE-001 distinction (§1),
superseding the earlier alias assumption now that this artifact is ingested.
</content>
