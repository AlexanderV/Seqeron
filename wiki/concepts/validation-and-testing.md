---
type: concept
title: "Validation and testing strategy"
tags: [testing, validation]
sources:
  - README.md
  - ALGORITHMS_CHECKLIST_V2.md
  - docs/ADVANCED_TESTING_CHECKLIST.md
  - docs/sonar-gate-plan.md
source_commit: 7bad3df2460ec31b7b1a3be1e605c11198342268
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:layered-architecture
      source: advanced-testing-checklist
      evidence: "Architecture Testing (ArchUnitNET) ... Enforces layer boundaries (Core -> Analysis -> IO) ... 22 rules across all 13 modules"
      confidence: high
      status: current
---

# Validation and testing strategy

Seqeron's development-time correctness story: 22,000+ executed test cases across 258 algorithm units (roughly 3.8× more test code than product code), green on .NET 10 with warnings-as-errors and CI-gated, plus a documented per-unit validation campaign. The "warnings-as-errors" half of that gate is its own discipline — see the [[build-quality-gate]] (SonarAnalyzer ratcheted to blocking under `TreatWarningsAsErrors`).

## Ten complementary methodologies

Each catches a different class of defect, with a per-algorithm checklist under `docs/checklists`. The [[advanced-testing-checklist]] rates all ten by applicability, current coverage, effort, and P0–P3 priority against the completed-unit set:

- **Property-based** (FsCheck) — invariant violations across generated inputs.
- **Metamorphic** — wrong outputs when the exact answer is unknown, via relations like `revcomp(revcomp(x)) == x`.
- **Fuzzing** — crashes/edge cases from malformed or adversarial input.
- **Mutation** (Stryker.NET) — weak tests, by seeding deliberate bugs.
- **Snapshot / approval** (Verify) — unintended changes to complex outputs.
- **Algebraic** — broken identity/inverse/idempotence/commutativity laws.
- **Architecture** (ArchUnitNET) — layering/dependency-rule drift (see [[layered-architecture]]).
- **Differential** — divergence from a reference implementation.
- **Combinatorial / pairwise** — interaction bugs across parameter spaces.
- **Characterization** — regressions during refactoring, by pinning behaviour.

Coverage across these ten is uneven and mostly aspirational. As of the [[advanced-testing-checklist]] baseline (2026-03-19), only **architecture testing** was marked complete; the two highest-priority techniques (property-based and metamorphic) were the biggest gap — property files existed for ~22 areas but many units lacked specific invariants, and 72 completed units had no metamorphic relations. Fuzzing and differential testing existed only for SuffixTree, not Genomics.

## Validation campaign

Beyond tests, the README describes a per-unit internal validation campaign: a findings register, a published limitations / operating-envelope document, literature-traced parameters, and one report per unit under `docs/Validation/reports`. This is the evidence base behind the runtime [[scientific-rigor]] guarantees — but it is internal and self-validated, with the caveats spelled out in [[research-grade-limitations]].

The campaign is tracked unit-by-unit in the [[test-unit-registry]] (364 units, 255 complete at the ingested revision), where each unit must clear the [[definition-of-done]]. Note the "258 algorithm units" figure from the README and the "364 test units" figure from [[algorithms-checklist-v2]] count different things — algorithm implementations versus tracked test units (255 done + 109 proposed) — and are not in conflict.
