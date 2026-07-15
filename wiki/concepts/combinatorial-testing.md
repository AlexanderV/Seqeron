---
type: concept
title: "Combinatorial / pairwise testing — covering parameter interactions minimally"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/09_COMBINATORIAL_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: combinatorial-testing-checklist
      evidence: "Priority P3. Generates a minimal set of test cases covering all parameter pairs (or t-tuples) — effective for algorithms with several configuration parameters where exhaustive enumeration is impractical (>100 combinations)."
      confidence: high
      status: current
---

# Combinatorial / pairwise testing

Combinatorial (pairwise) testing generates the **minimal set of test cases that covers every
pair** (or every t-tuple) of parameter values, on the empirical observation that most defects are
triggered by the interaction of *two* parameters, not four or five. For an algorithm with several
configuration knobs whose full cross-product runs into the hundreds, pairwise coverage collapses
the test count by an order of magnitude while still exercising every two-way interaction. This is
the lowest-priority (**P3**) member of the [[validation-and-testing]] program; the checklist
record is [[combinatorial-testing-checklist]].

## When it applies, and the tools

Applicability rule from the checklist: an algorithm takes **≥ 3 parameters**, the parameters have
discrete values or ranges, and the full enumeration exceeds **~100 combinations**. Complexity
banding: **Low = ≤ 2 params, Med = 3 params, High = ≥ 4 params**. Tooling: PICT (Microsoft),
AllPairs, and NUnit's `[Combinatorial]` / `[Pairwise]` attributes with `[Values]` / `[Range]`.

## Coverage

From a starting point of **zero**, the campaign reached **193 / 255 complete** with **65 marked ✗
not applicable** (≤ 2 parameters, so pairwise adds nothing) and an estimated **~900 pairwise test
cases** generated. Priority split: **15** high (≥ 4 params, > 50 full combinations), **52**
medium (3 params), **19** low. Combinatorial testing is orthogonal to the correctness
methodologies — it improves the *input coverage* of the existing property and example tests
rather than adding a new oracle, which is why it and [[characterization-testing]] sit at P3.
