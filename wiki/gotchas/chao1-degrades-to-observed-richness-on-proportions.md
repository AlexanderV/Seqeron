---
type: gotcha
title: "alpha_diversity Chao1 silently returns observed richness for non-integer (proportional) abundances"
tags: [metagenomics, gotcha]
mcp_tools:
  - alpha_diversity
sources:
  - docs/algorithms/Metagenomics/Alpha_Diversity.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# alpha_diversity Chao1 silently returns observed richness for non-integer (proportional) abundances

**The trap.** Chao1 estimates **unseen** richness from **singletons** (`f₁`) and **doubletons**
(`f₂`) — `Ŝ = S_obs + f₁²/(2·f₂)` — and singletons/doubletons only exist for **integer count
data**. The implementation therefore **gates Chao1 on the input looking like counts**: if the
abundances are non-integer / proportional (relative abundances, TPM, fractions), it **returns
`ObservedSpecies` for `Chao1Estimate`** instead of an extrapolated estimate. The output record shape
is unchanged, so **nothing signals the degradation**. Likewise, with no singletons (`f₁ = 0`),
Chao1 = `S_obs`.

**Why it bites.** Relative-abundance tables are the **common** metagenomics form. Feed one in and
`Chao1Estimate` silently equals `ObservedSpecies`: you believe you have an unseen-richness estimator
but you are reading back observed richness. Comparing "Chao1" across samples where some are counts
and some are proportions compares **two different quantities**.

**What to rely on instead.** Pass **raw integer counts** to get a real Chao1 extrapolation. If you
only have proportions, do **not** interpret `Chao1Estimate` as an estimate (it is observed
richness). `ShannonIndex`, `SimpsonIndex`, `InverseSimpson`, and `PielouEvenness` are defined on
proportions and are unaffected. Full model: [[alpha-diversity]]; cross-sample comparison uses
[[beta-diversity]].
