---
type: gotcha
title: "distance_matrix p-distance is uncorrected; JC69/K2P return +∞ at saturation"
tags: [phylogenetics, gotcha]
mcp_tools:
  - distance_matrix
  - build_tree_from_matrix
sources:
  - docs/algorithms/Phylogenetics/Distance_Matrix.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# distance_matrix p-distance is uncorrected; JC69/K2P return +∞ at saturation

**The trap.** The **p-distance** model is the raw Hamming **proportion** (`differences /
comparableSites`) with **no correction** for unseen (multiple / back) substitutions. It
systematically **underestimates** evolutionary divergence as sequences diverge. The corrected
models go the other way at biological **saturation** and become undefined: **JC69 → +∞ when
p ≥ 3/4**; **K2P → +∞ when a log argument goes non-positive**.

**Why it bites.**

- **Trees from p-distances are distorted.** Saturated, genuinely distant pairs are reported as
  closer than they are, so branch lengths and topology from `build_tree_from_matrix` are wrong for
  divergent input.
- **Switching to JC69/K2P can emit `+∞` / undefined entries** on very divergent sequences, which
  break downstream matrix math (NJ/UPGMA) if not handled as a saturation sentinel.
- Choosing the substitution model is a **scientific decision**, not a formatting toggle.

**What to rely on instead.** Use **JC69 or K2P** (or a richer model) for divergent sequences and
handle the `+∞` saturation value explicitly; reserve **p-distance** for closely related sequences
where the correction is negligible. Full model: [[evolutionary-distance-matrix]]; the tree builders
are [[tree-statistics]] / [[tree-comparison-metrics]]. Pipeline view:
[[comparative-genomics-pipeline-silent-traps]].
