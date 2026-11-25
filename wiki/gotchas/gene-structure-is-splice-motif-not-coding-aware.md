---
type: gotcha
title: "predict_gene_structure pairs GT-AG splice motifs — it is not coding/frame-aware, and splicing is intron removal"
tags: [splicing, gotcha]
mcp_tools:
  - predict_gene_structure
  - predict_introns
sources:
  - docs/algorithms/Splicing/Gene_Structure_Prediction.md
  - docs/Validation/review-2026-09/B11.md
source_commit: 8b85d3b5496273e86f3be2eef7f720a6ce7a257b
created: 2026-07-20
updated: 2026-11-25
---

# predict_gene_structure pairs GT-AG splice motifs — it is not coding/frame-aware, and splicing is intron removal

> **Updated 2026-11-25 (B11 F32–F36, F45).** The default `model=greedy` is still a splice-motif model (not coding-aware), but:
> exons + introns now **tile** the input (no exon is dropped; spliced sequence = input minus selected introns, F34); U12 typing needs
> both a U12 donor and a branch-site acceptor (F32); the intron score is a fixed (donor + acceptor + branch)/3 term (F33);
> `DetectAlternativeSplicing` uses SUPPA/rMATS SE / A5SS / A3SS / MXE definitions (F35). A **coding-aware** alternative now exists:
> `model=splice_dp` (`PredictGeneStructureDp`) — an exact MaxEntScan splice + open-reading-frame DP (no premature in-frame stop;
> returns the CDS) — but with no trained length/composition priors its precision is low (23–34 %). Details: [[gene-structure-prediction-intron-exon]];
> report `docs/Validation/review-2026-09/B11.md`.


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
