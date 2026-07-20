---
type: gotcha
title: "predict_chou_fasman returns Pα/Pβ/Pt propensity profiles, not H/E/C secondary-structure calls"
tags: [protein-features, gotcha]
mcp_tools:
  - predict_chou_fasman
sources:
  - docs/algorithms/Statistics/Secondary_Structure_Prediction.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# predict_chou_fasman returns Pα/Pβ/Pt propensity profiles, not H/E/C secondary-structure calls

**The trap.** `predict_chou_fasman` emits, for each window, the **arithmetic mean of the Chou-Fasman
Pα / Pβ / Pt propensities** — a **sliding-window propensity profile** (the same `N − W + 1`
window-mean machinery as the Kyte-Doolittle hydropathy profile). It does **not** run the full
Chou-Fasman **nucleation + extension** procedure and does **not** assign a per-residue
helix/strand/coil (**H/E/C**) state.

**Why it bites.** If you expect a secondary-structure **string** (like DSSP/PSIPRED output) you'll
instead get three propensity tracks. Deciding "this residue is a helix" by taking `argmax(Pα, Pβ)`
per position is **not** the Chou-Fasman algorithm — the real method requires seeded nucleation
windows and directional extension with break rules, which this profile does not perform.

**What to rely on instead.** Read the output as **conformational propensity tracks** for comparing
regions; for actual H/E/C state assignment use a real predictor (PSIPRED / JPred / a full
Chou-Fasman implementation). Full model: [[protein-secondary-structure-chou-fasman]]; the analogous
profile primitive is [[hydrophobicity-gravy-and-profile]].
