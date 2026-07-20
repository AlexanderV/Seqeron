---
type: gotcha
title: "taxonomic_profile is read-tally relative abundance over classified reads — not MetaPhlAn marker abundance"
tags: [metagenomics, gotcha]
mcp_tools:
  - taxonomic_profile
sources:
  - docs/algorithms/Metagenomics/Taxonomic_Profile.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# taxonomic_profile is read-tally relative abundance over classified reads — not MetaPhlAn marker abundance

**The trap.** `taxonomic_profile` aggregates the **per-read** classifications from
[[taxonomic-classification]] into **relative-abundance** estimates per taxon (kingdom / phylum /
genus / species). It is a **read-tally** profiler — abundance = share of classified reads — **not**
the clade-specific-marker model of MetaPhlAn, and the fractions are computed over the **classified**
reads only (unclassified reads don't create entries).

**Why it bites.** Read-tally abundance inherits read-count biases: taxa with **larger genomes or
higher rRNA/gene copy number** contribute more reads and read **higher** than their true cell
fraction, so the profile is not directly comparable to a MetaPhlAn marker-based profile. Because the
denominator is classified reads, a sample with many unclassified reads yields relative abundances
that **do not reflect the whole community** — they sum to 1 over what was classifiable.

**What to rely on instead.** Interpret the output as relative composition of the **classified**
fraction; report the classified/unclassified split, and use a marker-based profiler if you need
genome-size-normalized abundances. Downstream diversity summaries are [[alpha-diversity]] /
[[beta-diversity]]. Full model: [[taxonomic-profile]].
