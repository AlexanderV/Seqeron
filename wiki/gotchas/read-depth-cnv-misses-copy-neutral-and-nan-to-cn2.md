---
type: gotcha
title: "segment_copy_number is read-depth-only: blind to copy-neutral events, and NaN windows become neutral CN2"
tags: [structural-variants, gotcha]
mcp_tools:
  - segment_copy_number
sources:
  - docs/algorithms/StructuralVar/Copy_Number_Variation.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# segment_copy_number is read-depth-only: blind to copy-neutral events, and NaN windows become neutral CN2

**The trap.** This is a **read-depth** CNV caller (Yoon 2009): it assumes **read depth is
proportional to copy number** and segments the GC-adjusted mean depth. Because it only sees depth, it
is **blind to copy-neutral rearrangements** (balanced translocations, inversions, copy-neutral LOH),
and a **NaN / no-signal window is replaced with Neutral (CN 2)** — a silent default, not a missing
call. Status *Simplified*, single-sample, **not for clinical use**.

**Why it bites.** Depth-only calling inherits the classic biases: **mappability and GC** confound the
ratio (the optional GC correction only normalizes to the median, it doesn't fix low-mappability
regions), balanced events are invisible, and a stretch of uncovered/masked windows reads back as
**normal CN 2** rather than "unknown". Interpreting a CN-2 region as confidently diploid can hide a
copy-neutral event or a coverage gap.

**What to rely on instead.** Combine with breakpoint-level evidence for balanced/copy-neutral events
([[breakpoint-detection-split-reads]], [[discordant-pair-sv-detection]]); treat NaN→CN2 windows as
low-confidence, and account for GC/mappability. For tumor allele-specific copy number use the ASCAT
path ([[allele-specific-copy-number-ascat]]). Full model: [[read-depth-cnv-segmentation]].
