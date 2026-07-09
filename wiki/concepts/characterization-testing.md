---
type: concept
title: "Characterization testing — pinning current behavior before a refactor"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/10_CHARACTERIZATION_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: characterization-testing-checklist
      evidence: "Priority P3. Golden-master tests that pin the system's current 'as-is' behaviour before a refactor — 'не перевіряють коректність — перевіряють незмінність поведінки' (they don't check correctness, they check behavioural invariance)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:snapshot-testing
      source: characterization-testing-checklist
      evidence: "The Snapshots/ tests effectively play a similar role, but characterization tests are specific to refactoring — generated on-demand before a risky change and discarded after, whereas snapshots are a standing regression guard."
      confidence: high
      status: current
---

# Characterization testing

Characterization tests (golden-master tests) **pin the system's current behaviour "as-is"** so
that a refactor can be proven behaviour-preserving. Crucially, they **do not check correctness —
they check invariance**: they capture whatever the code does today, correct or not, and fail on
any deviation. They are generated **on demand before a risky change** and are throwaway, not a
standing suite. This is a **P3** member of the [[validation-and-testing]] program; the checklist
record is [[characterization-testing-checklist]].

## When and how

Apply before: replacing an algorithm with a new implementation; optimising (Span-based, SIMD);
extracting code into a separate module; or changing an API's parameters or return types. The
process: (1) generate a set of inputs (corner cases + typical cases), (2) record current outputs
as a golden master, (3) refactor, (4) run — **any divergence fails**, (5) review the diff — an
intentional change is approved, a regression is fixed.

## Coverage and relation to snapshot testing

Formal coverage is **0** — characterization tests are on-demand by nature, so there is no
standing count to complete. The checklist notes that the `Snapshots/` files ([[snapshot-testing]])
"effectively play a similar role," and the two are near-duplicates in mechanism (serialize output,
diff against a committed golden file). The distinction is *intent and lifecycle*: snapshots are a
permanent regression guard run every build; characterization tests are a temporary safety net
erected around one refactor. Both differ fundamentally from the correctness methodologies
([[property-based-testing]], [[metamorphic-testing]], [[algebraic-testing]]) in that they assert
**sameness, not correctness** — which is exactly why they carry the program's lowest priority.
