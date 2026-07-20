---
type: gotcha
title: "find_shared_motifs counts sequences containing a motif, not total occurrences (first occurrence only)"
tags: [motif, gotcha]
mcp_tools:
  - find_shared_motifs
sources:
  - docs/algorithms/Motif_Discovery/Shared_Motifs.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# find_shared_motifs counts sequences containing a motif, not total occurrences (first occurrence only)

**The trap.** Following the RSAT convention, `find_shared_motifs` takes **"only the first occurrence
of each sequence into account"**: a motif is credited to a sequence when it appears **at least
once**, and repeated occurrences **within the same sequence are not counted**.

**Why it bites.** If you read the output as a **total-occurrence count** or an abundance, you will be
wrong. A motif that occurs 10× in one sequence and 1× in another scores the **same** as one that
occurs 1× in each — both are "present in 2 sequences". Ranking motifs by how "shared" they are is a
**presence/absence count across input sequences**, not a frequency of hits.

**What to rely on instead.** Treat the result as **sequence-level presence/absence** — how many
input sequences contain the motif. For per-sequence or genome-wide **occurrence counts**, use an
occurrence-enumerating search such as `find_exact_motif` / known-motif search
([[known-motif-search]]). Full model: [[shared-motifs]].
