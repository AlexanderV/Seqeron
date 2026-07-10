---
type: source
title: "Evidence: RNA-STRUCT-001 (Secondary Structure Prediction — Nussinov base-pair maximization + constraints)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-STRUCT-001-Evidence.md
sources:
  - docs/Evidence/RNA-STRUCT-001-Evidence.md
source_commit: bb82b7ec80bbbf5750e53616ccc60df7b45c010c
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-STRUCT-001

The validation-evidence artifact for test unit **RNA-STRUCT-001** — **Secondary Structure Prediction**
(area `RnaStructure`): the **top-level structure-prediction umbrella** on `RnaSecondaryStructure` whose
canonical method is `Predict(sequence)` (a **Nussinov** base-pair-maximizing dynamic program), plus
**constraint folding** (`PredictWithConstraints`) and **dot-bracket I/O conversion** (`ToDotBracket` /
`FromDotBracket`). It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept is
[[rna-secondary-structure-prediction]]. [[test-unit-registry]] tracks the unit.

## Relationship to RNA-MFE-001 (naming correction)

This is a **distinct sibling** of [[rna-minimum-free-energy-folding|RNA-MFE-001]], **not** an alias.
Prior RNA ingests, before this artifact was ingested, assumed the generic id *RNA-STRUCT-001* denoted
the same folder as RNA-MFE-001. The now-ingested artifact shows otherwise: RNA-STRUCT-001's headline
algorithm is **Nussinov base-pair maximization** (`Predict`) — plus constraints and notation
conversion — whereas RNA-MFE-001 is the physical **Turner-2004 kcal/mol Zuker–Stiegler MFE folder**
(`CalculateMinimumFreeEnergy` / `PredictStructure`). Both units exercise the same `RnaSecondaryStructure`
class and share its Zuker MFE machinery: **deviation D5** of *this* artifact added
`CalculateMfeStructure` / `PredictStructureMfe` (a Zuker–Stiegler traceback over the same V/W/WM
matrices as the scalar MFE), and the greedy `Predict` / `PredictStructure` is retained as a fast path.
So the two ids are **two test units in one area**, exercising overlapping code, with different headline
algorithms — not two names for one unit.

## What this file records

- **Algorithm(s):**
  - **Nussinov & Jacobson (1980)** DP — the O(n³) time / O(n²) space dynamic program that **maximizes
    the number of base pairs**; cannot detect pseudoknots (correct by design).
  - **MFE DP (weighted-pair Nussinov variant):** Nussinov-style DP with **weighted pair scores** —
    Watson-Crick **−2.0**, wobble **−1.0** — maximizing weighted pair count, **not** thermodynamic MFE.
    Results indicate **relative stability, not physical energy units** (contrast RNA-MFE-001's kcal/mol).
  - **Stem-loop energy model:** Turner 2004 nearest-neighbor stacking (NNDB), terminal AU/GU penalty
    (+0.45 per helix end), hairpin-loop initiation (sizes 3–30, Jacobson–Stockmayer beyond 30), special
    tri/tetra/hexaloop total-energy overrides, first-mismatch bonuses (UU/GA −0.9, GG −0.8), all-C loop
    penalty (+1.5 for 3 nt, 0.3n+1.6 for n>3).
  - **Zuker-style DP traceback** (`CalculateMfeStructure` / `PredictStructureMfe`, deviation D5, added
    2026-06-23) — the reconstructed structure's energy equals `CalculateMinimumFreeEnergy` for the same
    input (asserted across hairpin, multi-stem, multiloop cases). Pseudoknotted optima remain excluded
    (the O(n³) recurrences are pseudoknot-free by construction).

- **Method coverage:** `Predict(sequence)` (canonical, complete prediction), `PredictWithConstraints(seq,
  constraints)` (forced base pairs — Mathews et al. 2004 constrained-DP class), `ToDotBracket(structure)`
  (structure → notation), `FromDotBracket(notation)` (notation → structure).

- **Authoritative sources:** Nussinov & Jacobson 1980 (*PNAS* 77(11):6309, original O(n³) base-pair-max
  DP); Zuker & Stiegler 1981 (*NAR* 9(1):133, MFE + W/V matrices + traceback); MIT 6.047 Lecture 08
  (Washietl 2012 — explicit Zuker F/C/M/M¹ recurrences, `F_ij = min{F_{i+1,j}, min_k C_ik + F_{k+1,j}}`,
  `C_ij = min{hairpin, interior, multiloop}`); Turner 2004 / NNDB (nearest-neighbor parameters);
  Mathews et al. 2004 (*PNAS* 101(19):7287, chemical-modification/constraint DP); Wikipedia (Nucleic
  acid structure prediction, Nussinov algorithm, Nucleic acid secondary structure). Secondary: Rosetta
  Code, ViennaRNA docs, MFOLD server.

- **Structural motifs classified:** stem (helix), hairpin loop, internal loop, bulge, multi-loop
  (junction), pseudoknot (crossing pairs — detection only, not predicted).

- **Reference oracles:**
  - Simple hairpin `GGGGAAAACCCC` → `((((....))))`, 4-bp stem + 4-nt loop, MFE negative.
  - GNRA tetraloop `GCGCGAAACGCGC` (G-GAAA-C closing) → GA first-mismatch bonus −0.9; with G-C closing
    the standard model + GA bonus applies (NNDB special table is C-G-closing only).
  - tRNA-like `GCGGAUUUAGCUCAGUUGG…GAAUUCGCA` (72 nt) → multiple stem-loops, cloverleaf, valid
    dot-bracket.
  - Poly-A `AAAAAAAAAAAA` → no base pairs, MFE = 0, all-dots dot-bracket (A cannot pair with A).

- **Base-pairing invariants (shared [[rna-base-pairing]] rule):** Watson-Crick {A-U, U-A, G-C, C-G} +
  wobble {G-U, U-G} pair; A-A/A-G/A-C/U-U/U-C/C-C/G-G do not.

- **Structural invariants:** dot-bracket balance (openers = closers); MFE sign ≤ 0 (Nussinov weighted
  pair score); WC stacking always negative; loop initiation always positive (destabilizing); selected
  stem-loops non-overlapping; each base in ≤ 1 pair. Probability invariants: structure P ∈ [0,1], MFE
  structure highest probability, Boltzmann `exp(−E/RT)/Z`.

- **Pseudoknot detection oracles:** non-crossing (0,5)+(1,4) → none; crossing (0,6)+(3,9) → detected;
  nested (0,10)+(2,8)+(4,6) → none (the `i<k<j<l` crossing test shared with [[rna-pseudoknot-detection]]).

- **Edge / corner cases:** empty `""` and null → empty structure, MFE = 0; too-short `"GC"` → no
  stem-loop; case-insensitive (`gggaaaaccc` ≡ uppercase); invalid characters handled gracefully or
  rejected; minimum hairpin = 3-bp stem + 3-nt loop (steric floor); wobble-only stems fold if wobble
  enabled; poly-A / poly-U → no structure.

## Deviations and assumptions

- **Resolved:** D1 single-nucleotide-bulge degeneracy `−RT·ln(states)` added (verified NNDB Example 1:
  3 C's → −0.68 kcal/mol); D2 dangling ends (model d2) added to the multiloop WM recurrence; D5 Zuker
  traceback `CalculateMfeStructure`/`PredictStructureMfe` (above).
- **Open (blocked):** D3 1×2 internal-loop `int21` lookup (2,304 entries) and D4 2×2 `int22` lookup
  (36,864 entries) — too large for inline static tables; a generic initiation + asymmetry + mismatch
  model is used instead. Would require an external data file.
- **Design decisions (not deviations):** no pseudoknot prediction (O(n³) DP cannot; correct by design);
  single sequence only (no comparative/covariance); minimum loop size 3 (NNDB steric constraint).

**No source contradictions** — Nussinov & Jacobson 1980, Zuker & Stiegler 1981, MIT 6.047, Turner 2004 /
NNDB, and Mathews 2004 agree on the DP decompositions, the base-pair-max vs weighted-pair objectives,
and the worked energies. The one **wiki reconciliation** is the RNA-STRUCT-001 ≠ RNA-MFE-001 distinction
noted above (the earlier alias assumption is superseded now that this artifact is ingested).
</content>
</invoke>
