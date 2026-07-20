---
type: gotcha
title: "find_discordant_pairs gives SV signatures and approximate intervals, not base-pair breakpoints"
tags: [structural-variants, gotcha]
mcp_tools:
  - find_discordant_pairs
sources:
  - docs/algorithms/StructuralVar/SV_Detection.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# find_discordant_pairs gives SV signatures and approximate intervals, not base-pair breakpoints

**The trap.** This is the **paired-end-mapping (PEM)** member of the SV family (Medvedev, Stanciu &
Brudno 2009): a **discordant read pair** (anomalous insert size or orientation) is a **signature** of
an SV within the insert-size envelope — it flags *that a rearrangement is near here and of this
orientation type*, **not the exact breakpoint**. Base-pair resolution comes from within-read
**split-read** evidence ([[breakpoint-detection-split-reads]]) or breakpoint assembly, not from the
discordant pair itself.

**Why it bites.** Reading a discordant-pair cluster as a precise breakpoint is wrong by up to the
library insert size. Orientation tells you the **SV type** (deletion / duplication / inversion /
translocation signature), but the coordinates are an **interval**, and PEM alone cannot resolve the
junction sequence or microhomology.

**What to rely on instead.** Cluster discordant pairs for SV **discovery and typing**, then refine
the breakpoint with split reads / breakpoint assembly for base-pair coordinates. Copy-number-changing
events also show up in [[read-depth-cnv-segmentation]]. Full model:
[[discordant-pair-sv-detection]]. Research-grade, [[research-grade-limitations]].
