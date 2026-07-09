---
type: concept
title: "Algorithm validation evidence artifacts"
tags: [validation, testing]
sources:
  - docs/Evidence/ALIGN-GLOBAL-001-Evidence.md
source_commit: 46d4efa2e08a672c942aa455eeb8b724705081e3
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-global-001-evidence
      evidence: "Evidence Artifact: ALIGN-GLOBAL-001 ... Online Sources ... Test Dataset ... Deviations and Assumptions ... References"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:definition-of-done
      source: align-global-001-evidence
      evidence: "Deviations and Assumptions: None ... the implementation follows the standard Needleman–Wunsch linear gap penalty model exactly as given in the Wikipedia pseudocode"
      confidence: medium
      status: current
---

# Algorithm validation evidence artifacts

Each algorithm [[test-unit-registry|test unit]] has a per-unit **Evidence artifact** under
`docs/Evidence/<UnitID>-Evidence.md`. These are the literature-traced source record behind
the [[definition-of-done]]'s "Evidence documented" criterion and the
[[validation-and-testing]] campaign: they pin exactly which external references and worked
examples the implementation and its tests are validated against.

## Templated structure

Every evidence file follows the same shape:

1. **Header** — Test Unit ID, algorithm name, date collected.
2. **Online sources** — Wikipedia / primary-literature URLs with access dates and the key
   extracted points (definitions, recurrences, complexity, worked examples).
3. **Test dataset** — the canonical worked example(s) with exact parameters and expected
   outputs, used as the oracle / differential test fixture.
4. **Deviations and assumptions** — where the implementation departs from (or exactly
   follows) the reference, plus API-contract behaviours (null/empty handling) that sit
   outside the algorithm spec.
5. **References** — primary literature and encyclopedic citations.

Because these files are near-templated across the ~213 documented units, the wiki keeps
**one** shared page for the pattern (this page) plus a concise per-file source summary
(e.g. [[align-global-001-evidence]]). An individual algorithm gets its own concept page
only when it is itself distinct and wiki-worthy — for example
[[global-alignment-needleman-wunsch]].
