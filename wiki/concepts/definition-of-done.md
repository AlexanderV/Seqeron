---
type: concept
title: "Definition of Done (test units)"
tags: [validation, testing]
sources:
  - ALGORITHMS_CHECKLIST_V2.md
source_commit: 6a14170477c9472c0be07e3b7c7f7123e31eddcf
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: algorithms-checklist-v2
      evidence: "Definition of Done (DoD) ... Required Criteria: TestSpec created, Tests written, Branch coverage >= 80%, Edge cases covered, Tests pass, Evidence documented"
      confidence: high
      status: current
---

# Definition of Done (test units)

The acceptance bar a [[test-unit-registry|test unit]] must clear before it is marked complete (`☑`), defined in [[algorithms-checklist-v2]]. It is what makes "done" objective rather than a judgement call.

## Required criteria

1. **TestSpec created** — `TestSpecs/<TestUnitID>.md`.
2. **Tests written** — `*.Tests/<Class>_<Method>_Tests.cs`.
3. **Branch coverage ≥ 80%** — per the coverage report.
4. **Edge cases covered** — null, empty, boundary, error.
5. **Tests pass** — CI green.
6. **Evidence documented** — PR/commit link in the registry.

## Quality and complexity criteria

Beyond the six, tests must be independent (order-independent), deterministic (no unseeded randomness), named `Method_Scenario_ExpectedResult`, structured Arrange-Act-Assert, with one assert per logical check. Algorithms that are **O(n²) or higher** additionally require a **property-based test for the invariant** and a **recorded performance baseline**.

This DoD operationalizes the project's [[scientific-rigor|rigor]] at development time and is the gate that populates the [[validation-and-testing]] numbers. Note that some units ship ahead of full validation (see [[research-grade-limitations]]), so a `☑` reflects this bar while `☐` states can mean "code present, DoD/re-validation not yet met."
