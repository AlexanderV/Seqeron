---
type: concept
title: "Architecture testing (ArchUnitNET) — preventing dependency drift at the IL level"
tags: [testing, validation, methodology, architecture]
sources:
  - docs/checklists/07_ARCHITECTURE_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: architecture-testing-checklist
      evidence: "Priority P2, framework ArchUnitNET. 22 module-level rules (all ☑) checking layer boundaries, no cyclic dependencies, no System.IO in Core, immutable Result/DTO types, and placement/naming invariants — verified at the IL level."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:layered-architecture
      source: architecture-testing-checklist
      evidence: "The 22 rules enforce the Core/IO layer boundaries across all 13 modules and forbid cyclic dependencies — the executable guard on the up-only dependency layering described by layered-architecture."
      confidence: high
      status: current
---

# Architecture testing (ArchUnitNET)

Architecture testing prevents **architectural drift** by asserting **structural rules over the
compiled IL** — which module may depend on which, naming conventions, type placement — as
executable tests that fail the build when violated. Unlike the other nine methodologies its rules
apply to **modules (projects), not individual algorithms**, so a new algorithm is *automatically*
covered by its module's rules without any new test. This is the executable enforcement of the
[[layered-architecture]] Dependency Rule. It is a **P2** member of the
[[validation-and-testing]] program and, per the [[advanced-testing-checklist]], was the **first
of the ten to reach complete**; the checklist record is [[architecture-testing-checklist]].

## The 22 rules

`Architecture/ArchitectureTests.cs` holds **22 rules, all complete**, over **13 modules**:

- **Layer boundaries** — `Core` and `IO` boundaries to all 13 modules; the up-only dependency
  direction is enforced, not merely documented.
- **No cyclic dependencies** between modules.
- **No `System.IO` in `Core`** — I/O concerns stay out of the domain core.
- **Immutable `Result`/DTO types** — result and data-transfer types cannot expose mutable state.
- **Placement / naming invariants** — parsers live only in `IO`; `Core` holds no algorithm
  classes; namespace matches assembly.

The module layout the rules police: `Core` (base types, sequences), `Analysis` (89 algorithms —
k-mer, repeat, motif, disorder, GC-skew, complexity, comparative), `Alignment` (16 — NW/SW/MSA +
assembly), `IO` (9 — parsers + Phred quality), `Annotation` (37 — ORF/gene/promoter, miRNA,
splicing, epigenetics, variants, SV, transcriptome), `MolTools` (23 — CRISPR/primer/probe/
restriction/codon), plus `Phylogenetics`, `Population`, `Chromosome`, and `Metagenomics`.

## Why it scales differently

Because the rules bind modules rather than units, architecture testing's **22/22 completeness is
stable** as the library grows — the opposite of the per-unit methodologies whose coverage must
chase every new algorithm. It is the structural counterpart to the correctness-oriented siblings
([[property-based-testing]], [[metamorphic-testing]], [[mutation-testing]]) and pairs with the
[[build-quality-gate]] (warnings-as-errors, Sonar-blocking) as the two always-green guardrails.
