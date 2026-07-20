---
type: gotcha
title: "predict_coiled_coils is a heptad a/d-occupancy heuristic, not a COILS/Marcoil probability"
tags: [protein-features, gotcha]
mcp_tools:
  - predict_coiled_coils
sources:
  - docs/algorithms/ProteinMotif/Coiled_Coil_Prediction.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# predict_coiled_coils is a heptad a/d-occupancy heuristic, not a COILS/Marcoil probability

**The trap.** `PredictCoiledCoils` is a **heuristic, sequence-only** predictor: for each sliding
window it computes the **fraction of heptad a/d core positions filled by a hydrophobic residue**
(`{I, L, V}`) — a register-occupancy score. It is **not** the COILS position-specific scoring matrix
or a Marcoil HMM probability.

**Why it bites.** The score is an occupancy fraction, not a calibrated coiled-coil probability, and
it keys on a single heptad frame's hydrophobic-core occupancy — so it can miss coiled-coils with
imperfect registers and flag any hydrophobic-periodic stretch. Comparing the value to a COILS
`P(cc)` or ranking candidates as if it were one is over-reading.

**What to rely on instead.** Use it as a fast periodicity screen; confirm hits with COILS / Marcoil /
DeepCoil for decision-grade calls. It is one of several sequence-only windowed protein heuristics
(alongside [[intrinsic-disorder-prediction-top-idp]] and Chou-Fasman propensity) — each scores a
different signal, none is a trained model. Full model: [[coiled-coil-prediction]].
