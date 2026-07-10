---
type: source
title: "Evidence: PROTMOTIF-CC-001 (Coiled-coil prediction — heptad a/d hydrophobic-core occupancy)"
tags: [validation, protein]
doc_path: docs/Evidence/PROTMOTIF-CC-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-CC-001-Evidence.md
source_commit: 0003535e3857fc7a8321c06da24972ffd5a14383
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-CC-001

The validation-evidence artifact for test unit **PROTMOTIF-CC-001** — **coiled-coil
prediction** (`ProteinMotifFinder.PredictCoiledCoils`): score each sliding window by the
fraction of heptad **a/d core positions** occupied by a hydrophobic residue `∈ {I, L, V}`,
maximised over the **seven heptad registers**, and report contiguous high-scoring spans as
coiled-coil regions. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the model, contract, invariants
and worked oracle are synthesized in [[coiled-coil-prediction]]. See [[test-unit-registry]]
for how units are tracked.

## What this file records

- **Online sources:**
  - **Mason JM & Arndt KM 2004**, *Coiled coil domains* (rank 1, ChemBioChem 5(2):170–176) —
    heptad `[abcdefg]n` with hydrophobic a/d and polar/charged e/g; a/d burial ("almost
    complete van der Waals contact") drives helix association; interactions among a, d, e, g
    give most structural specificity. Supplies the **≥3-heptad multi-heptad** requirement.
  - **Wikipedia "Coiled coil"** (rank 4, citing Mason & Arndt) — fixes the hydrophobic-core
    residue set verbatim: a and d "often being occupied by **isoleucine, leucine, or valine**"
    → **{I, L, V}**.
  - **Wikipedia "Heptad repeat"** (rank 4, citing Chambers 1990) — `HPPHCPC` pattern, a/d
    hydrophobic; leucine zippers have predominantly **L at d** (confirms L is canonical).
  - **Lupas, Van Dyke & Stock 1991**, *Predicting coiled coils from protein sequences*
    (rank 1, Science 252:1162–1164) — the COILS method: **7 heptad registers** (frames) and a
    canonical **28-residue window** (4 heptads); "196 preliminary scores per residue"
    (7 frames × 28 window positions), score = max over windows/registers. The full COILS
    **21×20 PSSM was not retrievable**, so it is deliberately **NOT** implemented; only the
    fully-specified register + 28-window machinery is reused.
- **Datasets (documented oracles):** a/d occupancy is a **closed-form count**, not a tabulated
  empirical constant — `score = (#a/d positions with residue ∈ {I,L,V}) / (#a/d positions)`,
  max over 7 registers:
  - `(LAALAAA)` repeated (L at both a and d every heptad) → all a/d = L → **1.0**.
  - All-glycine → no a/d ∈ {I,L,V} → **0.0**.
  - `(LAAAAAA)` repeated (L at a only) → half of a/d filled → **0.5** (threshold boundary).
- **Corner cases / failure modes:** sub-window sequence cannot be scored → no prediction;
  register unknown a priori so all 7 frames must be tried and the best taken; single isolated
  hydrophobic positions are not coiled coils (multi-heptad requirement); no {I,L,V} → no
  regions; off-frame coiled coil still found via the 7-register max; lowercase recognised;
  null/empty → empty.

## Deviations and assumptions

**Deviations:** the COILS position-specific frequency matrix is **omitted** (weights not
retrievable from authoritative sources; users needing PSSM P-scores should use COILS/Paircoil2)
— the a/d occupancy fraction is used instead, so e/g electrostatic-edge residues do not
contribute and scores are occupancy fractions in `[0,1]`, not COILS P-scores. Two
**assumptions**, both source-anchored parameter choices:

1. **Hydrophobic-core set = {I, L, V}** — exactly the set named verbatim by the source; A/M/F
   coiled coils may score lower, accepted to avoid untraceable constants.
2. **Defaults window = 28, MinRegion = 21, threshold = 0.5** — window 28 and 7 registers from
   Lupas 1991; MinRegion 21 (3 heptads) from Mason & Arndt's multi-heptad requirement;
   threshold 0.5 = "predominantly hydrophobic" (> half of a/d filled). All caller-overridable.

Recommended coverage (MUST): perfect heptad (L at every a/d) → one region, score 1.0 spanning
the sequence; no {I,L,V} → no regions; sequence shorter than window → empty; off-frame coil
found via 7-register max; region shorter than MinRegion (21) rejected. SHOULD: half-occupancy
(L at a only) → exactly 0.5 (boundary); null/empty → empty. COULD: case-insensitivity. No
contradictions among sources.
