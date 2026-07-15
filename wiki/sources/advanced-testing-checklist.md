---
type: source
title: "Advanced Testing Checklist — technique effectiveness analysis"
tags: [validation, testing]
doc_path: docs/ADVANCED_TESTING_CHECKLIST.md
source_commit: 7bad3df2460ec31b7b1a3be1e605c11198342268
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Advanced Testing Checklist — technique effectiveness analysis

A point-in-time (dated **2026-03-19**) audit of how well ten advanced testing techniques fit Seqeron.Genomics, scoring each on applicability, current coverage, effort, and priority, then recording the concrete gap and a per-technique checklist. It is the "which techniques, and where are we short" companion to the [[validation-and-testing]] strategy and its [[algorithms-checklist-v2]] registry.

## What it adds over the strategy page

[[validation-and-testing]] names the ten methodologies; this document *rates and gap-analyses* them against the then-completed unit set. Each technique gets a 1–5★ applicability score, an existing-coverage count, an effort estimate, and a P0–P3 priority, plus its current state and the specific units still uncovered.

| Technique | Applicability | Coverage at ingest | Priority |
|---|:---:|---|:---:|
| Property-based (FsCheck) | ★★★★★ | 22 property files | P0 |
| Metamorphic | ★★★★★ | 7 units (`MetamorphicTests.cs`, 18+ MRs) | P0 |
| Mutation (Stryker) | ★★★★☆ | 2 source files (MotifFinder, RepeatFinder) | P1 |
| Algebraic | ★★★★☆ | 0 dedicated (some implicit in property files) | P1 |
| Snapshot / approval (Verify) | ★★★★☆ | ~20 snapshot files | P1 |
| Architecture (ArchUnitNET) | ★★★☆☆ | ☑ 22 rules, all 13 modules | P2 |
| Differential | ★★★☆☆ | 0 (SuffixTree safe-vs-unsafe only) | P2 |
| Fuzzing | ★★★☆☆ | 0 (SuffixTree header-corruption only) | P2 |
| Characterization | ★★☆☆☆ | 0 explicit (snapshots cover it) | P3 |
| Combinatorial / pairwise | ★★☆☆☆ | 0 explicit | P3 |

The two P0 techniques (property-based, metamorphic) are the headline gap: property files exist for ~22 areas but many units within them lack specific invariant tests, and 72 of the completed units have no metamorphic relations despite strong natural ones (alignment symmetry, GC preservation under complement, diversity monotonicity). See `docs/ADVANCED_TESTING_CHECKLIST.md` for the full per-technique rationale, and the ten per-technique checklists it links under `docs/checklists/01_…10_*.md`.

## Notable specifics

- **Architecture testing is the one technique marked done** — `Architecture/ArchitectureTests.cs` carries 22 IL-level rules enforcing the [[layered-architecture]] boundaries (Core → Analysis → IO), no inter-module cycles, no `System.IO` in Core, immutable Result/DTO types, and placement/naming invariants. Documented in `docs/checklists/07_ARCHITECTURE_TESTING.md`.
- **Fuzzing and differential testing both exist only for SuffixTree** (a separate library), never for Genomics; parsers (FASTA/FASTQ/BED/VCF/GFF/GenBank/EMBL) and sequence validation are the named prime fuzz targets.
- **Combinatorial testing** is aimed at the multi-knob algorithms — primer design, codon optimization, CRISPR, restriction digests — where ≥3 independent parameters make full enumeration expensive.
- An **execution roadmap** (Part 2) orders the work P0→P3: fill property gaps + 7 new files, ~15 new metamorphic MR classes, per-module Stryker runs, algebraic laws folded into property files, ~31 missing snapshots, fuzz files for 7 parsers, differential pairs, then on-demand characterization and pairwise suites.

## Caveats and contradictions

- **Snapshot in time, superseded counts.** The scope line says "**79 completed test units**," but Part 2 says the checklists cover "**ALL 86 completed algorithms**" — an internal 79-vs-86 discrepancy. Both are lower than the **255 completed** recorded later in [[algorithms-checklist-v2]] (version 2.5, 2026-02-12 registry snapshot). The numbers reflect different scopes and dates; treat this document's coverage figures as a 2026-03-19 baseline, not current state.
- **Aspirational, not attained.** Except for architecture testing (☑), the priorities and roadmap describe *planned* coverage. The gaps it lists reinforce [[research-grade-limitations]]: the testing program is deep but still has large uncovered surface in the P0 techniques.

## Where this fits

- [[validation-and-testing]] — this document is the effectiveness/gap analysis behind that strategy's ten methodologies.
- [[algorithms-checklist-v2]] — the unit registry whose completed units this audit scores.
- [[layered-architecture]] — enforced by the one technique already complete (architecture testing).
- [[research-grade-limitations]] — the outstanding coverage gaps are part of the research-grade caveat.
