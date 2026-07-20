---
type: gotcha
title: "PredictGenes emits every ORF as a CDS — no overlap resolution, no gene ranking"
tags: [annotation, gotcha]
mcp_tools:
  - predict_genes
sources:
  - docs/algorithms/Annotation/Gene_Prediction.md
source_commit: ec9209f6a267e376a0cd93f5b2e02d3576035966
created: 2026-07-11
updated: 2026-07-20
---

# PredictGenes emits every ORF as a CDS — no overlap resolution, no gene ranking

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
