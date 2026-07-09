---
type: source
title: "Checklist 09: Combinatorial / Pairwise Testing"
tags: [validation, testing, methodology]
doc_path: docs/checklists/09_COMBINATORIAL_TESTING.md
sources:
  - docs/checklists/09_COMBINATORIAL_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 09: Combinatorial / Pairwise Testing

The **P3** per-unit checklist for combinatorial/pairwise testing — a 255-row table tracking
minimal parameter-interaction coverage. Synthesized in the concept [[combinatorial-testing]];
part of the [[validation-and-testing]] program tracked in the [[test-unit-registry]].

## What this file records

- **Purpose:** generate a minimal set of cases covering all parameter pairs (or t-tuples);
  effective when full enumeration is impractical (> 100 combinations).
- **Applicability rule:** ≥ 3 parameters, discrete values/ranges, > 100 full combinations.
  Complexity banding Low = ≤ 2 params, Med = 3, High = ≥ 4.
- **Tools:** PICT, AllPairs, NUnit `[Combinatorial]` / `[Pairwise]` / `[Values]` / `[Range]`.
- **Starting point:** 0 pairwise tests. **Per-unit summary:** **193 ☑, 65 ✗ not applicable**
  (≤ 2 params), ~900 pairwise cases estimated; priority split 15 high (≥ 4 params), 52 medium,
  19 low.

## Deviations and contradictions

None. Orthogonal to the correctness methodologies — it improves *input coverage* of existing
tests, not the oracle, hence P3 alongside [[characterization-testing]].
