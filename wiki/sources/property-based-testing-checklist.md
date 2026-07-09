---
type: source
title: "Checklist 01: Property-Based Testing (FsCheck)"
tags: [validation, testing, methodology]
doc_path: docs/checklists/01_PROPERTY_BASED_TESTING.md
sources:
  - docs/checklists/01_PROPERTY_BASED_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 01: Property-Based Testing (FsCheck)

The **P0** per-unit checklist for property-based testing — a 258-row status table mapping every
test unit to the invariants it must satisfy and the `Properties/*.cs` file that verifies them.
Synthesized in the concept [[property-based-testing]]; part of the ten-methodology
[[validation-and-testing]] program tracked in the [[test-unit-registry]].

## What this file records

- **Framework:** FsCheck + FsCheck.NUnit; discipline that every genomic algorithm expresses ≥ 1
  invariant, checked over hundreds of generated inputs.
- **Invariant legend:** R (Range), S (Symmetry), I (Idempotence/Involution), M (Monotonicity),
  P (Preservation), RT (Round-trip), D (Determinism).
- **Per-unit table:** all **258 units ☑ complete**; each row lists the specific invariants
  (e.g. `SEQ-GC-001`: GC% ∈ [0,100], complement preserves GC%) and the property file.
- **Summary:** 258 complete, 0 not started; 4 new property files added (Chromosome,
  Epigenetics, Oncology), ~15 existing files extended.

## Deviations and contradictions

None internal. Against the older [[advanced-testing-checklist]] (2026-03-19 baseline, which
called property-based testing the biggest gap with "property files existed for ~22 areas but many
units lacked specific invariants"), this checklist records the **completed** end-state — a
temporal progression, not a contradiction.
