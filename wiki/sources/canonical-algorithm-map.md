---
type: source
title: "Canonical Algorithm Map (docs/algorithms canonical IDs, folder buckets, legacy baselines)"
tags: [meta, coverage, governance]
doc_path: docs/algorithms/CANONICAL_MAP.md
sources:
  - docs/algorithms/CANONICAL_MAP.md
source_commit: f776181cd75c204d4932cac5e47d04285d03052f
ingested: 2026-07-13
created: 2026-07-13
updated: 2026-07-13
---

# Canonical Algorithm Map

A **meta / index** document under `docs/algorithms/` — not a per-algorithm spec. It is the
project's authority for **algorithm identity**: it collapses duplicate test-unit IDs onto one
canonical ID, normalizes the `docs/algorithms/**` folder buckets, and lists the legacy/baseline
methods that are kept on purpose. Its job is to prevent duplicate units and split ownership so that
one concept has exactly one canonical home. It governs the same `docs/algorithms/**` tree that the
wiki's coverage ledger — [[backlog]] — reconciles against concept pages, and its canonical IDs are
the identities registered in the test-unit registry [[algorithms-checklist-v2]] and validated in
the [[algorithm-validation-evidence]] artifacts. Reference doc: `docs/algorithms/CANONICAL_MAP.md`.

## What it governs

- **Canonical ID aliases.** Three alias IDs are consolidated onto a canonical ID (same fixture /
  behavior): `SEQ-COMPOSITION-001` → `SEQ-STATS-001`; `SEQ-TM-001` → `SEQ-THERMO-001`;
  `GENOMIC-TANDEM-001` → `REP-TANDEM-001` (unified under repeats). Aliases stay searchable but must
  point to the canonical ID.
- **Taxonomy / folder-name aliases.** Variant folder names map to one canonical bucket:
  `Molecular_Tools` → `MolTools`; `PopGen` → `Population_Genetics` (short prefix kept in IDs only);
  `K-mer_Analysis` → `K-mer`; `RNA_Secondary_Structure` / `RNA_Structure` → `RnaStructure`.
- **Canonicalization rules.** One canonical ID per concept; aliases resolve to it; exactly one
  canonical behavior document per concept (domain pages link, never copy); legacy-baseline methods
  are retained and explicitly labeled `legacy-baseline`.
- **Legacy / baseline methods.** Kept as standard comparators, not defaults, each with a short
  baseline note linking back: **UPGMA** (`Phylogenetics/Tree_Construction.md`; molecular-clock
  assumption — prefer Neighbor-Joining), **Jukes-Cantor / Kimura-2P**
  (`Phylogenetics/Distance_Matrix.md`), **Chi-square Hardy-Weinberg**
  (`Population_Genetics/Hardy_Weinberg_Test.md`; exact tests preferable for sparse counts),
  **Nussinov-style RNA baseline** (`RnaStructure/Minimum_Free_Energy.md`; simplified energy model),
  **OLC assembly** (`Assembly/Overlap_Layout_Consensus.md`; often less practical than de Bruijn).

## Status recorded by the map

Folder normalization for the four canonical buckets above is **complete** — the alias folders
(`K-mer_Analysis`, `Molecular_Tools`, `PopGen`, `RNA_Structure`, `RNA_Secondary_Structure`) have
been merged and every inbound reference repointed. Remaining consolidation candidates named but not
yet done: `Motif_Analysis`, `Sequence_Comparison`, `Genomic_Analysis`, `Extended_Annotation`,
`Extended_Assembly`.

## Relation to the wiki's coverage infrastructure

- [[backlog]] tracks *coverage* (which `docs/algorithms/**` doc is synthesized by which concept
  page); this map tracks *identity* (which ID/folder is canonical). They are complementary: the
  backlog already flags `CANONICAL_MAP.md` and `docs/algorithms/README.md` as index/map docs, not
  algorithm units. This source page satisfies that navigational need without creating per-algorithm
  concepts.
- The canonical IDs here are the unit IDs registered in [[algorithms-checklist-v2]] and carried
  through the validation pipeline ([[validation-ledger]], [[algorithm-validation-evidence]]).
