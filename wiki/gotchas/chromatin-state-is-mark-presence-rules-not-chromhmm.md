---
type: gotcha
title: "predict_chromatin_state is a presence/absence mark-combination rule set, not ChromHMM on read counts"
tags: [epigenetics, gotcha]
mcp_tools:
  - predict_chromatin_state
  - annotate_histone_modifications
  - find_accessible_regions
sources:
  - docs/algorithms/Epigenetics/Chromatin_State_Prediction.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# predict_chromatin_state is a presence/absence mark-combination rule set, not ChromHMM on read counts

**The trap.** Chromatin state here is a function of the **set of present histone marks**, **not their
magnitudes**. It assigns a state from which marks are flagged present (rule-based combinatorics —
e.g. H3K4me3 + H3K27me3 → bivalent/poised TSS), and it does **not** perform ChromHMM's
**Poisson-background binarization from raw read counts** or a learned emission/transition model.

**Why it bites.** You must supply **binary present/absent marks**; feeding signal magnitudes and
expecting the tool to threshold them itself is wrong — magnitude is ignored once a mark is present.
And the state labels come from a fixed rule table, not a trained ChromHMM/Segway model, so genome
segmentations will not match a ChromHMM run and won't discover novel states beyond the coded rules.

**What to rely on instead.** Binarize marks upstream (your own peak-calling / thresholding), then use
this for the **state-assignment logic**; for a trained genome segmentation from read counts run
ChromHMM/Segway. Full model: [[chromatin-state-prediction]]; sequence-only CpG signal is
[[cpg-island-detection]]. Research-grade, [[research-grade-limitations]].
