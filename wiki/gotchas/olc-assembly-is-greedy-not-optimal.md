---
type: gotcha
title: "assemble_olc is greedy (order-dependent, error-free-read) — not an exact shortest-common-superstring"
tags: [assembly, gotcha]
mcp_tools:
  - assemble_olc
  - find_overlap
  - find_all_overlaps
sources:
  - docs/algorithms/Assembly/Overlap_Layout_Consensus.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# assemble_olc is greedy (order-dependent, error-free-read) — not an exact shortest-common-superstring

**The trap.** The exact shortest-common-superstring is **NP-complete** (Compeau 2011), so
`assemble_olc` uses a **greedy maximal-overlap layout** (greedy-SCS) plus concatenation consensus.
Greedy layout is **suboptimal** — the result can **depend on merge order** (the same read set can
give different superstring lengths) and is a **heuristic, not the optimal** assembly. The published
numeric oracles assume **error-free reads** with **exact-match (identity 1.0)** overlaps.

**Why it bites.** Expecting a deterministic-optimal or reference-grade assembly from `assemble_olc`
is wrong: a different input ordering can change the contig set, and real reads with sequencing
errors won't overlap under the exact-match path, so error-containing data assembles poorly or not at
all. This is a *correctness reference* for the OLC paradigm, not a production assembler with error
modelling.

**What to rely on instead.** Treat `assemble_olc` as a validated OLC demonstrator on clean reads;
error-correct reads first ([[kmer-spectrum-error-correction]]) and don't rely on a single greedy run
being globally optimal. `find_overlap` / `find_all_overlaps` expose the overlap graph if you want to
inspect the layout. Full model: [[overlap-layout-consensus-assembly]]; the pairwise merge step is
[[contig-merge-overlap-collapse]]; QC contiguity via [[assembly-statistics]]. Pipeline view (this is
the root of the assembly-error cascade): [[comparative-genomics-pipeline-silent-traps]].
