---
type: concept
title: "Strictly-layered library architecture"
tags: [architecture]
sources:
  - README.md
source_commit: 260735c56cf01fd76968956e28281f7678fff716
created: 2026-07-09
updated: 2026-07-17
---

# Strictly-layered library architecture

Seqeron is organized as a strictly-layered set of packages where dependencies only ever point *up* the levels — never sideways within a level, never downward — a rule enforced by architecture tests (ArchUnitNET). In the README's dependency graph an arrow `A --> B` means **B depends on A**.

## The levels

- **Substrate** — `SuffixTree` (Ukkonen + persistent), the fast substring-query engine.
- **Level 0** — `Infrastructure`.
- **Level 1** — `Core`.
- **Level 2** — `IO`, `Alignment`, `Analysis`.
- **Level 3** — `Annotation`, `Phylogenetics`, `Population`, `Metagenomics`, `MolTools`, `Chromosome`, `Oncology`.
- **Level 4** — `Reports`.
- **Meta-package** — `Seqeron.Genomics`, which aggregates all modules and is the single `using` namespace.

Every module also references `Core` + `Infrastructure`; those universal edges are elided in the README graph for readability.

## Why it matters

The dependency rule keeps the ~250 algorithms cohesive under one namespace and makes layering drift a test failure rather than a latent design rot. It is the structural counterpart to the runtime discipline described in [[scientific-rigor]], and the module boundaries here are the natural targets for future `module`/`api` wiki pages ingested from `docs/**`.

The repository layout that realizes these levels (`src/Seqeron/Algorithms`, `src/Seqeron/Mcp`, etc.) is detailed in `README.md`.
