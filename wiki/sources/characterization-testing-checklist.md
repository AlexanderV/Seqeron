---
type: source
title: "Checklist 10: Characterization Testing"
tags: [validation, testing, methodology]
doc_path: docs/checklists/10_CHARACTERIZATION_TESTING.md
sources:
  - docs/checklists/10_CHARACTERIZATION_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 10: Characterization Testing

The **P3** checklist for characterization (golden-master) testing — pinning current behaviour
before a refactor. Synthesized in the concept [[characterization-testing]]; part of the
[[validation-and-testing]] program.

## What this file records

- **Purpose:** capture the system's current "as-is" behaviour before a refactor — checks
  behavioural **invariance, not correctness**; used on-demand.
- **When to apply:** before replacing an algorithm, optimising (Span/SIMD), extracting a module,
  or changing an API's parameters/return types.
- **Process:** generate inputs (corner + typical) → record outputs as golden master → refactor →
  run (any divergence fails) → review diff (intentional → approve, regression → fix).
- **Coverage:** **0** standing count (on-demand by nature); the `Snapshots/` files
  ([[snapshot-testing]]) effectively play a similar role, but characterization tests are
  refactor-specific and throwaway.

## Deviations and contradictions

None. Near-duplicate of [[snapshot-testing]] in mechanism (golden master); the distinction is
lifecycle — permanent regression guard vs temporary refactor safety net.
