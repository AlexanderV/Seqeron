---
type: source
title: "Checklist 08: Differential Testing"
tags: [validation, testing, methodology]
doc_path: docs/checklists/08_DIFFERENTIAL_TESTING.md
sources:
  - docs/checklists/08_DIFFERENTIAL_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 08: Differential Testing

The **P2** per-unit checklist for differential testing — a 255-row table tracking cross-checks
between two independent implementations. Synthesized in the concept [[differential-testing]];
part of the [[validation-and-testing]] program tracked in the [[test-unit-registry]].

## What this file records

- **Purpose:** compare outputs of two independent implementations on identical inputs to surface
  subtle implementation bugs; the second implementation is the oracle.
- **Strategy legend:** ALT (alternative algorithm), BRUTE (brute-force reference for small
  inputs), REF (hand-computed / reference library), DUAL (two in-project implementations).
- **Starting point (recorded):** `SafeVsUnsafeDifferentialTests` for SuffixTree, plus
  **PROTMOTIF-HMM-001** — Plan7 local+glocal DP vs an independent path-enumeration brute force
  (`Plan7ProfileHmm_ForwardBackwardDifferential_Tests`).
- **Per-unit table + summary:** **107 / 255 ☑, 151 not started**; ranked ~25 high-value
  (ALT/BRUTE feasible), ~35 medium (REF), ~26 lower (needs DUAL re-implementation).

## Deviations and contradictions

None. Partial by design — differential testing needs a *second* correct implementation, the most
expensive precondition among the ten, so it is deliberately incomplete like [[snapshot-testing]].
