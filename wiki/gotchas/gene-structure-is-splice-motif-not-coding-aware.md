---
type: gotcha
title: "predict_gene_structure pairs GT-AG splice motifs — it is not coding/frame-aware, and splicing is intron removal"
tags: [splicing, gotcha]
mcp_tools:
  - predict_gene_structure
  - predict_introns
sources:
  - docs/algorithms/Splicing/Gene_Structure_Prediction.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# predict_gene_structure pairs GT-AG splice motifs — it is not coding/frame-aware, and splicing is intron removal

**The trap.** This is a **splice-motif** model (status *Simplified*): it finds donor/acceptor sites
and pairs them into introns under the **GT-AG rule** (plus the `GC-AG` and `AT-AC`/U12 variants),
using a position-level non-overlap check. It does **not** analyse coding potential, reading frame,
or ORFs — the exon/intron architecture is derived purely from splice-site motifs and scores. And
`GenerateSplicedSequence` is defined as **intron removal**, *not* as concatenating only coding exons.

**Why it bites.** Treating the output as a trained gene model (GeneMark/AUGUSTUS-style, coding-aware)
is wrong: a GT…AG pair inside a UTR or a non-coding region is still paired as an intron, and a
"single-exon gene" simply means no qualifying `GT-AG`/`AT-AC` pair was found — not a coding-potential
verdict. The spliced sequence you get back is the transcript with introns removed, which is **not**
the CDS.

**What to rely on instead.** For coding-potential / frame evidence combine with
[[coding-potential-hexamer-score]] and ORF finding ([[open-reading-frame-detection]]); use the
individual splice-site units for "is *this* a real splice site". Tune `minIntronLength` /
`maxIntronLength` / `minScore` to your organism. Full model:
[[gene-structure-prediction-intron-exon]].
