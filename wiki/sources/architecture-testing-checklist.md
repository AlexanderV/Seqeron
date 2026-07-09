---
type: source
title: "Checklist 07: Architecture Testing (ArchUnitNET)"
tags: [validation, testing, methodology, architecture]
doc_path: docs/checklists/07_ARCHITECTURE_TESTING.md
sources:
  - docs/checklists/07_ARCHITECTURE_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 07: Architecture Testing (ArchUnitNET)

The **P2** checklist for architecture testing — 22 module-level rules (not per-algorithm) over 13
modules. Synthesized in the concept [[architecture-testing]]; part of the
[[validation-and-testing]] program and the executable guard on [[layered-architecture]].

## What this file records

- **Framework:** ArchUnitNET; rules verified at the **IL level**, applied to modules (projects),
  so new algorithms are auto-covered by their module's rules.
- **The 22 rules (all ☑):** Core/IO layer boundaries to all 13 modules; no cyclic module
  dependencies; no `System.IO` in Core; immutable Result/DTO types; placement/naming invariants
  (parsers only in IO, Core without algorithm classes, namespace = assembly).
- **Module map + algorithm counts:** Core 9, Analysis 89, Alignment 16, IO 9, Annotation 37,
  MolTools 23, plus Phylogenetics/Population/Chromosome/Metagenomics.
- **Summary:** 22 rules complete, 0 remaining, 13 modules covered. Dated 2026-06-27.

## Deviations and contradictions

None. Consistent with the [[advanced-testing-checklist]] note that architecture testing was the
**first** of the ten methodologies to reach complete — module rules don't scale with algorithm
count, so completeness is stable.
