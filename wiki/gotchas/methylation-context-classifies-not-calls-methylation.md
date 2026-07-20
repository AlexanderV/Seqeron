---
type: gotcha
title: "find_methylation_sites classifies cytosine context (CpG/CHG/CHH) from sequence — it does not call methylation"
tags: [epigenetics, gotcha]
mcp_tools:
  - find_methylation_sites
sources:
  - docs/algorithms/Epigenetics/Methylation_Analysis.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# find_methylation_sites classifies cytosine context (CpG/CHG/CHH) from sequence — it does not call methylation

**The trap.** This is a **sequence-only** classifier: it partitions each cytosine into its
methylation **sequence context** — **CpG**, **CHG**, or **CHH** — from the DNA sequence alone. It
reports *where a cytosine could be methylated and in what context*, **not whether it is methylated**.
Calling actual methylation status requires bisulfite reads (that is
[[bisulfite-methylation-calling]], a different unit).

**Why it bites.** The name reads like a methylation caller. Feed it a genome and you get every
cytosine's context, not a methylome — treating the output as methylated positions would be wrong.
Also, non-CpG contexts (**CHG/CHH**) are biologically relevant mainly in **plants and stem cells**
(Lister 2009); in typical mammalian somatic DNA the meaningful context is CpG, so a large CHG/CHH
list is expected sequence bookkeeping, not evidence of non-CpG methylation.

**What to rely on instead.** Use this to annotate context and to drive the
`GenerateMethylationProfile` aggregator; to get β-values / methylation status call
[[bisulfite-methylation-calling]] on bisulfite reads. CpG *density* (island detection) is
[[cpg-island-detection]]. Full model: [[methylation-context-classification]].
