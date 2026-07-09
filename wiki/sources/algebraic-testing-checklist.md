---
type: source
title: "Checklist 06: Algebraic Testing"
tags: [validation, testing, methodology]
doc_path: docs/checklists/06_ALGEBRAIC_TESTING.md
sources:
  - docs/checklists/06_ALGEBRAIC_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 06: Algebraic Testing

The **P1** per-unit checklist for algebraic testing — a 258-row table mapping units to the
algebraic laws they must satisfy. Synthesized in the concept [[algebraic-testing]]; part of the
[[validation-and-testing]] program tracked in the [[test-unit-registry]].

## What this file records

- **Purpose:** verify explicit algebraic laws — a more formal, targeted form of
  [[property-based-testing]].
- **Law legend:** ID (Identity), COMM (Commutativity), ASSOC (Associativity), INV (Involution),
  IDEMP (Idempotence), RT (Round-trip/Isomorphism), DIST (Distributivity/Conservation), TRI
  (Triangle inequality).
- **Per-unit table:** **89 ☑ complete, 169 ✗ not applicable**, ~172 laws verified (≈2 per
  applicable algorithm).
- **Starting point:** laws existed implicitly in property files (involution in
  `SequenceProperties`, round-trip in `FastaRoundTripProperties`) but were not systematised.

## Deviations and contradictions

None. The large **not-applicable** count is intentional honesty — not every operation has an
algebraic law, so the methodology is partial by design rather than force-fit.
