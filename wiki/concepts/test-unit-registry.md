---
type: concept
title: "Test-unit registry and ID scheme"
tags: [validation, testing]
sources:
  - ALGORITHMS_CHECKLIST_V2.md
source_commit: 6a14170477c9472c0be07e3b7c7f7123e31eddcf
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: algorithms-checklist-v2
      evidence: "Total Test Units 364 ... Completed 255 ... each unit maps to a TestSpec and test file(s) in the Processing Registry"
      confidence: high
      status: current
---

# Test-unit registry and ID scheme

The unit of account for Seqeron's validation effort: every algorithm is a **test unit** with a stable, area-prefixed ID (e.g. `SEQ-GC-001`, `ALIGN-GLOBAL-001`, `CRISPR-OFF-001`, `ONCO-SIG-004`). The registry — defined in [[algorithms-checklist-v2]] — is how completion, evidence, and test artifacts are tracked per unit.

## The ID scheme

IDs are `<AREA>-<NAME>-<NNN>`. Areas partition the library: Composition, Matching, Repeats, MolTools, Annotation, K-mer, Alignment, Phylogenetic, PopGen, Chromosome, Metagenomics, Codon, Translation, FileIO, RnaStructure, MiRNA, Splicing, ProteinPred/ProteinMotif, Epigenetics, Variants, StructuralVar, Assembly, Transcriptome, ComparativeGenomics, PanGenome, and Oncology.

## Per-unit record

Each unit carries:

- a **Processing Registry** row — area, method count, **evidence** (Wikipedia + primary literature + reference tools), a `TestSpecs/<ID>.md`, and the test file(s);
- a **spec block** — canonical method, complexity, invariant, and enumerated edge cases;
- a **status** — `☑` completed, `☐` not started (with a distinct `☐`-pending-re-validation state for campaign-added units).

At the ingested revision: **364 units, 255 completed, 109 proposed**, and 44/57 classes covered. Every unit must clear the [[definition-of-done]]. This registry is the concrete mechanism behind the [[validation-and-testing]] strategy, and future ingests of individual `TestSpecs/` documents should link back here by unit ID. The literature-traced source record for each unit lives in a per-unit [[algorithm-validation-evidence|evidence artifact]] under `docs/Evidence/`.
