---
type: concept
title: "Coiled-coil prediction (heptad a/d hydrophobic-core occupancy)"
tags: [analysis, algorithm]
mcp_tools:
  - predict_coiled_coils
sources:
  - docs/Evidence/PROTMOTIF-CC-001-Evidence.md
  - docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md
source_commit: 335716dcdb5101ac27141e185a3e963cbc9fdb25
created: 2026-07-10
updated: 2026-07-16
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: protmotif-cc-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-CC-001 ... Algorithm: Coiled-Coil Prediction (heptad-repeat a/d hydrophobic-core detection)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:protein-low-complexity-seg
      source: protmotif-cc-001-evidence
      evidence: "Sibling sequence-only protein-feature heuristics in Seqeron's docs/algorithms protein groups; both scan a window over the amino-acid sequence and emit contiguous regions, but coiled-coil scores a/d hydrophobic periodicity while SEG scores Shannon-entropy low complexity — distinct algorithms"
      confidence: medium
      status: current
---

# Coiled-coil prediction (heptad a/d hydrophobic-core occupancy)

Predicting **coiled-coil** regions — bundles of α-helices wound around one another via a
seven-residue **heptad repeat** whose hydrophobic-core positions drive helix association. Seqeron's
`ProteinMotifFinder.PredictCoiledCoils` is a **heuristic, sequence-only** predictor: for each
sliding window it computes the fraction of heptad **a/d core positions** filled by a hydrophobic
residue `∈ {I, L, V}`, **maximised over the seven possible heptad registers**, and reports
contiguous high-scoring spans. Validated under test unit **PROTMOTIF-CC-001**; the validation
record is [[protmotif-cc-001-evidence]] and [[test-unit-registry]] tracks the unit. See
[[algorithm-validation-evidence]] for the artifact pattern.

This is the **first ingested unit of the protein-motif / ProteinMotif family** (coiled-coil,
signal-peptide, transmembrane-helix, domain prediction). It is a **distinct algorithm** from the
protein-disorder / features family — [[protein-low-complexity-seg|SEG low-complexity]] (Shannon
entropy over 20 residues), [[intrinsic-disorder-prediction-top-idp|intrinsic-disorder prediction]]
(TOP-IDP propensity), and [[morf-prediction-dip-in-disorder|MoRF prediction]] — all of which are
sequence-only windowed heuristics but score different signals. Coiled-coil detection uniquely keys
off the **a/d hydrophobic periodicity** of the α-helical heptad, not compositional bias or disorder
propensity.

## Heptad model and the a/d occupancy score

Coiled coils recur as `(abcdefg)n`. Folded into an α-helix (3.6 residues/turn), positions **a** and
**d** fall on the same face and form a hydrophobic stripe whose burial ("knobs into holes") against
another helix stabilises the assembly (Mason & Arndt 2004). For an uppercased sequence `S`, a window
of length `W` at index `i`, and register `r ∈ {0,…,6}`:

- Heptad position of residue index `k` in register `r`: `p(k,r) = (k − r) mod 7`.
- **Core positions** are `a = 0` and `d = 3`.
- A residue is a **hydrophobic-core residue** iff it is in **{I, L, V}** — the set named verbatim by
  the source ("a and d … often being occupied by isoleucine, leucine, or valine").
- Register occupancy: `occ(i,r) = (#core positions in window with residue ∈ {I,L,V}) / (#core positions)`.
- Window score: `score(i) = max over r of occ(i,r)` — the register is unknown a priori, so all seven
  frames are tried and the best is taken (Lupas 1991).

A window is coiled-coil if `score(i) ≥ threshold`. Contiguous coiled windows `[i₀ … i₁]` map to the
residue region `[i₀, i₁ + W − 1]`, reported (with `Score` = peak `score` in the run) only if its
length is at least the **minimum 3 heptads = 21 residues**.

## Parameters and defaults

| Parameter | Default | Role / source |
|-----------|---------|---------------|
| `windowSize` | **28** | sliding-window length = 4 heptads (Lupas 1991); sequences shorter than this yield nothing |
| `threshold` | **0.5** | minimum a/d hydrophobic-core occupancy = "predominantly hydrophobic" (> half of a/d positions) |
| min region | **21** | 3-heptad multi-heptad requirement (Mason & Arndt); shorter scoring runs are rejected |
| registers | **7** | the seven heptad frames tried per window (Lupas: "196 preliminary scores" = 7 frames × 28 positions) |

All are caller-overridable. `windowSize`, `threshold`, min-region and the 7 registers are the same
machinery COILS (Lupas 1991) uses; only the **COILS 21×20 position-specific scoring matrix is
deliberately omitted** — its weights were not retrievable from authoritative sources, so the fully
specified a/d occupancy fraction is used instead (see deviations).

## Invariants

- **INV-01** — every `Score ∈ [0, 1]` (a count ratio hydrophobic/core).
- **INV-02** — each region spans ≥ 21 residues (3 heptads).
- **INV-03** — `0 ≤ Start ≤ End ≤ n − 1`; regions are non-overlapping and increasing in `Start`.
- **INV-04** — sequences shorter than `windowSize` produce no regions (no full window exists).
- **INV-05** — a region exists only if some covering window scores ≥ threshold (max over 7 registers).

## Canonical oracle and corner cases

a/d occupancy is a **closed-form count**, not a tabulated empirical constant, so exact oracles are
derivable:

- `(LAALAAA)` repeated — L at both a and d every heptad → occupancy **1.0** in register 0. The
  5-heptad `"LAALAAA"×5` (35 aa) → one region `(Start=0, End=34, Score=1.0)`.
- All-glycine (or any sequence with no {I,L,V}) → occupancy **0.0** < threshold → **no regions**.
- `(LAAAAAA)` repeated — L at a only → half the a/d positions filled → **0.5** (threshold boundary).

Other cases: sequence shorter than `windowSize` → empty; an **off-frame** coiled coil is still found
via the 7-register max (INV-05); lowercase input recognised (uppercased first); null/empty → empty; a
scoring run shorter than the 21-residue minimum is rejected (INV-02). Complexity is **O(n·W·7)** time
(O(n) with the fixed constants W=28, 7 registers) / O(n) space.

Implementation shape (`ProteinMotifFinder.cs`): `PredictCoiledCoils(sequence, windowSize, threshold)`
computes the whole score profile once into an array, then scans it in a **single forward pass**, yielding
regions **lazily**. The per-window register max is a private `BestHeptadOccupancy(...)` and the run→region
map with the min-length filter is a private `BuildRegion(...)`. Because the algorithm is a per-position
**numeric window scan** and performs no exact-substring search, the repository suffix tree is **not
applicable** here — there is no occurrence enumeration, only windowed scoring.

## Deviations, assumptions, and scope

Two documented items. **Deviation:** the COILS PSSM is **not implemented** — scores are a/d occupancy
fractions in `[0,1]`, not COILS P-scores, and e/g electrostatic-edge residues do not contribute; users
needing PSSM-based probabilities should use the dedicated COILS / Paircoil2 tools. **Assumption:** the
hydrophobic-core set is limited to **{I, L, V}** — exactly the set named in the source, so
A/M/F-rich coiled coils may score lower; accepted to avoid fabricating untraceable constants. The
default window 28, min-region 21 and threshold 0.5 are source-anchored parameter choices, all
caller-overridable. A sequence-only heuristic: it detects the a/d hydrophobic periodicity but does not
model helix propensity, e/g electrostatics, oligomerization state, or register breaks/stutters, and may
over-predict generic hydrophobic-periodic sequences. A [[research-grade-limitations|research-grade]]
implementation.

## References

Mason J.M. & Arndt K.M. (2004) *ChemBioChem* 5(2):170–176 (heptad a/d burial, multi-heptad
requirement); Lupas A., Van Dyke M. & Stock J. (1991) *Science* 252(5009):1162–1164 (COILS: 7
registers, 28-residue window); Chambers P. et al. (1990) *J Gen Virol* 71(12):3075–3080 (heptad
`HPPHCPC`); Wikipedia "Coiled coil" / "Heptad repeat" ({I,L,V} core set, leucine-zipper L at d). Full
citations in [[protmotif-cc-001-evidence]] (do not duplicate here).

**Sharp edge:** [[coiled-coil-score-is-heptad-occupancy-heuristic]] — the score is a **heptad a/d-occupancy heuristic**, not a COILS/Marcoil probability.
