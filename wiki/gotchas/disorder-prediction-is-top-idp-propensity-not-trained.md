---
type: gotcha
title: "predict_disorder is a TOP-IDP propensity window, not a trained disorder predictor (IUPred/DISOPRED)"
tags: [protein-features, gotcha]
mcp_tools:
  - predict_disorder
  - disorder_propensity
  - is_disorder_promoting
sources:
  - docs/algorithms/ProteinPred/Disorder_Prediction.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# predict_disorder is a TOP-IDP propensity window, not a trained disorder predictor (IUPred/DISOPRED)

**The trap.** `predict_disorder` scores each residue as the **sliding-window average of the
normalized TOP-IDP propensity scale** (Campen et al. 2008). It is a **sequence-only propensity
heuristic**, not a trained/energy-based predictor like IUPred, DISOPRED, or PONDR — there is no
learned model and no calibrated probability behind the number.

**Why it bites.** Reading the per-residue value as a **probability of disorder** over-interprets a
composition propensity: the profile tracks disorder-promoting amino-acid content, so a compositionally
biased but folded region can score high, and the absolute cutoff between "ordered" and "disordered"
is a threshold on a propensity, not a calibrated decision boundary. Results will not match IUPred/
DISOPRED on the same sequence.

**What to rely on instead.** Use it as a fast propensity screen and compare **relative** regions of
one protein; for decision-grade IDR calls run a trained predictor. The raw primitives
(`GetDisorderPropensity` + Dunker classification) are separately validated if you want the
per-residue scale directly. Full model: [[intrinsic-disorder-prediction-top-idp]]. Research-grade,
[[research-grade-limitations]].
