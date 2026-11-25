---
type: gotcha
title: "PredictGenes emits every ORF as a CDS — no overlap resolution, no gene ranking"
tags: [annotation, gotcha]
mcp_tools:
  - predict_genes
sources:
  - docs/algorithms/Annotation/Gene_Prediction.md
  - docs/Validation/review-2026-09/B11.md
source_commit: 8b85d3b5496273e86f3be2eef7f720a6ce7a257b
created: 2026-07-11
updated: 2026-11-25
---

# PredictGenes emits every ORF as a CDS — no overlap resolution, no gene ranking

> **Updated 2026-11-25 (B11 F6, F44).** `PredictGenes` now emits **one gene per stop codon** (longest ORF = most upstream ATG/GTG/TTG; getorf-confirmed), not every nested start. Opt-in `startSelection` ProdigalScore / ProdigalDynamicProgramming (Prodigal V2.6.3 model trained on the input genome) resolves start choice and overlaps; overlap resolution/ranking in the default mode is still absent.


**The trap.** `GenomeAnnotator.PredictGenes(...)` looks like a gene finder, but it is
ORF-first *annotation scaffolding*, not a trained pipeline. It emits **every** qualifying ORF
(≥ `minOrfLength` amino acids, both strands) as its own annotation — including **overlapping
and nested candidates that share a stop codon**. There is no best-model selection, no overlap
suppression, and no ranking or filtering by promoter / RBS / coding-potential evidence
(`PredictGenes` never calls the RBS helper). Every record is labelled `Type = "CDS"` and
`Product = "hypothetical protein"`.

**Why it bites.** A caller expecting one gene per locus (GLIMMER / GeneMark-style output) will
instead get a redundant, inflated list — the same coding region can appear several times as
different in-frame start choices. Deduplication, overlap resolution, and any strength ranking
are the **caller's** job. Two further sharp edges compound it:

- **`minOrfLength` is in amino acids** (default 100), unlike the Analysis-layer
  `GenomicAnalyzer.FindOpenReadingFrames` whose `minLength` is in **nucleotides** — a unit trap
  when porting a threshold between the two ORF entry points.
- **RBS records are ORF-driven and forward-only** in the legacy `FindRibosomeBindingSites`, so a
  predicted reverse-strand gene may carry no ribosome-binding-site hit from that method.

**What to rely on instead.** For promoter evidence use the separate [[promoter-detection]]
helper; for coding-potential scoring use [[coding-potential-hexamer-score]]. There is no
in-repo codon-bias / organism-trained / intron-aware gene finder — the limitation is declared,
not a defect. Full model: [[prokaryotic-gene-prediction-rbs]] (test unit ANNOT-GENE-001,
primary spec `docs/algorithms/Annotation/Gene_Prediction.md`).
