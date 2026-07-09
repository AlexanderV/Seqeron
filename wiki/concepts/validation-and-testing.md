---
type: concept
title: "Validation and testing strategy"
tags: [testing, validation]
sources:
  - README.md
  - ALGORITHMS_CHECKLIST_V2.md
source_commit: 6a14170477c9472c0be07e3b7c7f7123e31eddcf
created: 2026-07-09
updated: 2026-07-09
---

# Validation and testing strategy

Seqeron's development-time correctness story: 22,000+ executed test cases across 258 algorithm units (roughly 3.8× more test code than product code), green on .NET 10 with warnings-as-errors and CI-gated, plus a documented per-unit validation campaign.

## Ten complementary methodologies

Each catches a different class of defect, with a per-algorithm checklist under `docs/checklists`:

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

## Validation campaign

Beyond tests, the README describes a per-unit internal validation campaign: a findings register, a published limitations / operating-envelope document, literature-traced parameters, and one report per unit under `docs/Validation/reports`. This is the evidence base behind the runtime [[scientific-rigor]] guarantees — but it is internal and self-validated, with the caveats spelled out in [[research-grade-limitations]].

The campaign is tracked unit-by-unit in the [[test-unit-registry]] (364 units, 255 complete at the ingested revision), where each unit must clear the [[definition-of-done]]. Note the "258 algorithm units" figure from the README and the "364 test units" figure from [[algorithms-checklist-v2]] count different things — algorithm implementations versus tracked test units (255 done + 109 proposed) — and are not in conflict.
