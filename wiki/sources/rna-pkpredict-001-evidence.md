---
type: source
title: "Evidence: RNA-PKPREDICT-001 (Pseudoknot structure prediction — canonical H-type, pknotsRG class)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-PKPREDICT-001-Evidence.md
sources:
  - docs/Evidence/RNA-PKPREDICT-001-Evidence.md
source_commit: b5cb44721acd570224bd0d138372f30e6f95b2a5
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-PKPREDICT-001

The validation-evidence artifact for test unit **RNA-PKPREDICT-001** — **pseudoknot structure
prediction** of the canonical **H-type** (pknotsRG class): predict the energetically optimal RNA
secondary structure that may contain a single simple recursive pseudoknot (two crossing helices).
It is one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the synthesizing concept is [[rna-pseudoknot-prediction]]. [[test-unit-registry]] tracks the
unit.

This is the **crossing-helix layer** of the RNA secondary-structure family — it goes **beyond** the
nested (pseudoknot-free) structure space of [[rna-minimum-free-energy-folding|MFE folding]] and
[[rna-partition-function-mccaskill|the McCaskill ensemble]], both of which enumerate only non-crossing
folds. A pseudoknot is exactly a pair of **crossing** helices that no nested predictor can represent.

## What this file records

- **Algorithm:** canonical **simple recursive pseudoknot** prediction — "two crossing helices with
  three intervening loops" (Reeder & Giegerich 2004, pknotsRG). The grammar motif is
  `knot = a ~~~ u ~~~ b ~~~ v ~~~ a' ~~~ w ~~~ b'`: segment *a* pairs with *a'*, *b* with *b'*, and
  *u*/*v*/*w* are the three loops separating the helices (each can fold internally). Full pknotsRG is
  **O(n⁴) time, O(n²) space**.

- **Authoritative sources:**
  - **Reeder & Giegerich (2004)** *BMC Bioinformatics* 5:104 (rank 1, DOI 10.1186/1471-2105-5-104,
    PMC514697) — the pknotsRG algorithm: canonical simple recursive pseudoknots, O(n⁴)/O(n²), the
    canonization rules, and the pseudoknot penalties.
  - **pknotsRG reference source** `Energy.lhs` / `Foldingspace.lhs` (rank 3, github.com/jensreeder/
    pknotsRG) — verbatim penalties and the confirmation that helices use the **same** Turner NN model
    as nested structures with **no** extra per-base-pair penalty inside the knot.
  - **Wikipedia — Pseudoknot** (rank 4, cites Rivas & Eddy 1999) — the H-type 5'→3' order and
    two-layer `()`/`[]` annotation.
  - **RCSB PDB 437D / Su et al. (1999)** *Nat Struct Biol* 6(3):285 (rank 1–5) — the BWYV
    frameshifting pseudoknot (28-nt real H-type) whose tertiary stabilization sits **outside** the NN
    secondary-structure model.
  - **Antczak et al. (2018)** *Bioinformatics* 34(8):1304 — the crossing condition `i < k < j < l`.

- **Energy model (sourced, not invented):** Turner nearest-neighbour stacking for **both** helices,
  plus dangling/coaxial terms, plus pseudoknot-specific penalties verbatim from `Energy.lhs`:
  - **pseudoknot initiation** (creating a new pseudoknot): **9.0 kcal/mol** (Reeder & Giegerich
    "setting this value to 9 kcal/mole performs better").
  - **unpaired nucleotide inside a pseudoknot loop:** **0.3 kcal/mol** each.
  - **base pair inside the pseudoknot:** **0.0** (no extra per-pair penalty; scored with the same NN
    model as nested pairs).

- **Canonization rules (define the class):**
  1. Both strands of a helix have equal length (`|a| = |a'|`, `|b| = |b'|`); helices have **no
     bulges**.
  2. Facing helices have **maximal extent** — extended to maximal Watson-Crick/GU length, so helix
     length is determined by the sequence, not searched independently.
  3. If two maximal helices would overlap, their boundary is fixed at an arbitrary point.

- **H-type geometry (Wikipedia / Rivas & Eddy):** 5'→3' order stem1-5' → loop1 → stem2-5' → loop2 →
  stem1-3' → loop3 → stem2-3'; half of one stem is intercalated between the two halves of the other,
  producing the crossing arrangement written as a **two-layer** dot-bracket (`()` + `[]`).

## Oracles and datasets

- **Designed canonical H-type (fully derivable):** `GGGGAACCCCAACCCCAAGGGG` (22 nt) — layout
  S1a=`GGGG`[0–3] · L1=`AA` · S2a=`CCCC`[6–9] · L2=`AA` · S1b=`CCCC`[12–15] · L3=`AA` ·
  S2b=`GGGG`[18–21]. Stem 1 = (0,15)(1,14)(2,13)(3,12) four G·C; stem 2 = (6,21)(7,20)(8,19)(9,18)
  four C·G. Crossing check: S2a (6–9) lies inside S1's span (0–15), S2b (18–21) lies outside →
  `i<k<j<l`. Two-layer dot-bracket **`((((..[[[[..))))..]]]]`**; pseudoknot ΔG strictly below the
  plain-MFE ΔG of the same sequence (knot accepted).
- **Plain hairpin (no spurious pseudoknot):** `GGGGAAAACCCC` → `HasPseudoknot == false`; returned
  structure and ΔG equal the plain MFE `((((....))))`.
- **BWYV real knot (documented NN-thermodynamic non-recovery):** `GGCGCGGCACCGUCCGCGGAACAAACGG`
  (28 nt, PDB 437D). The NN optimum is the pseudoknot-free stem-1 hairpin; the crystallographic knot
  is tertiary-stabilized (minor-groove triplex, ion coordination) and is **not** the MFE structure —
  an expected limit of all NN-only pseudoknot predictors, prevents over-claiming.

## Invariants and edge cases

- **MFE fallback bound:** for any sequence, `FreeEnergy ≤ CalculateMfeStructure(seq).FreeEnergy` — the
  predictor never returns a structure worse than the plain MFE (the always-available baseline).
- **No spurious knot:** the 9 kcal/mol initiation penalty exists so a pseudoknot is reported only when
  the two crossing helices' stabilization outweighs both the penalty and the best pseudoknot-free
  alternative. When a nested structure is more stable, `HasPseudoknot == false`.
- **Validity when a knot is returned:** every index in range, each position paired at most once, and
  `DetectPseudoknots` finds ≥1 genuine crossing (`i<k<j<l`, Antczak 2018; canonization rule 1).
- **Null / empty / too-short → empty pseudoknot-free structure** (no pairs, all dots, ΔG 0) — contract
  parity with `CalculateMfeStructure`.
- **DNA input** (T read as U) folds identically to the RNA spelling (module parity, COULD-test).

## Deviations and assumptions

1. **PARTIAL coverage of the pknotsRG class (documented, not an invented parameter).** The
   implementation realizes the canonical **single** H-type pseudoknot (two crossing helices + three
   internally-foldable loops) with the sourced energy model and penalties. The full pknotsRG O(n⁴)
   grammar additionally composes **recursively-nested** pseudoknots and **over-arching / multiple**
   knots within one structure; these are **NOT** implemented. Loops *u*/*v*/*w* fold with the existing
   pseudoknot-free MFE (`CalculateMinimumFreeEnergy`), consistent with "loops can fold internally" but
   without re-searching a second knot inside a loop.

**No source contradictions** — Reeder & Giegerich 2004, the pknotsRG `Energy.lhs` reference source,
the Wikipedia/Rivas & Eddy H-type geometry, and the BWYV structural record agree on the canonical
class, the 9.0 / 0.3 / 0.0 penalties, the same-NN-model helix scoring, and the documented
tertiary-stabilization limit. The only recorded item is the PARTIAL-coverage scope note.
