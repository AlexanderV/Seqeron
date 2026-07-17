---
type: concept
title: "Validation and testing strategy"
tags: [testing, validation]
sources:
  - README.md
  - ALGORITHMS_CHECKLIST_V2.md
  - docs/ADVANCED_TESTING_CHECKLIST.md
  - docs/sonar-gate-plan.md
  - docs/Validation/FINDINGS_REGISTER.md
  - docs/Validation/LIMITATIONS.md
  - docs/Validation/VALIDATION_LEDGER.md
source_commit: 260735c56cf01fd76968956e28281f7678fff716
created: 2026-07-09
updated: 2026-07-17
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

- **[[property-based-testing|Property-based]]** (FsCheck) — invariant violations across generated inputs. **P0.**
- **[[metamorphic-testing|Metamorphic]]** — wrong outputs when the exact answer is unknown, via relations like `revcomp(revcomp(x)) == x`. **P0.**
- **[[fuzzing|Fuzzing]]** — crashes/edge cases from malformed or adversarial input. **P2.**
- **[[mutation-testing|Mutation]]** (Stryker.NET) — weak tests, by seeding deliberate bugs. **P1.**
- **[[snapshot-testing|Snapshot / approval]]** (Verify) — unintended changes to complex outputs. **P1.**
- **[[algebraic-testing|Algebraic]]** — broken identity/inverse/idempotence/commutativity laws. **P1.**
- **[[architecture-testing|Architecture]]** (ArchUnitNET) — layering/dependency-rule drift (see [[layered-architecture]]). **P2.**
- **[[differential-testing|Differential]]** — divergence from a reference implementation. **P2.**
- **[[combinatorial-testing|Combinatorial / pairwise]]** — interaction bugs across parameter spaces. **P3.**
- **[[characterization-testing|Characterization]]** — regressions during refactoring, by pinning behaviour. **P3.**

Coverage was uneven at the [[advanced-testing-checklist]] baseline (2026-03-19) — only **architecture testing** was complete, and the two P0 techniques (property-based, metamorphic) were the biggest gap. The **per-methodology checklists under `docs/checklists`** record a much-advanced later state: [[property-based-testing]], [[metamorphic-testing]], and [[fuzzing]] are now **258/258**, [[architecture-testing]] 22/22, [[combinatorial-testing]] 193/255, [[mutation-testing]] all-files-≥-80% (2026-06-30), [[algebraic-testing]] 89 (+169 not-applicable), [[differential-testing]] 107/255. The real remaining gap is **[[snapshot-testing]] (37/255)** and the on-demand [[characterization-testing]] — so the program is deep but still uneven, not uniformly done.

## Validation campaign

Beyond tests, the README describes a per-unit internal validation campaign: a findings register, a published [[limitations|limitations / operating-envelope document]] (enforced at runtime — see [[operating-envelope-and-limitation-policy]]), literature-traced parameters, and one report per unit under `docs/Validation/reports`. This is the evidence base behind the runtime [[scientific-rigor]] guarantees — but it is internal and self-validated, with the caveats spelled out in [[research-grade-limitations]]. Each unit is worked through the [[validation-protocol]]: a **fresh session per unit** in a context deliberately separate from the implementer's, doing **Stage A (validate the description against external primary sources) before Stage B (validate the code realises it)**, and ending in exactly one of two states — ✅ CLEAN or 🔧 LIMITED.

The [[findings-register]] is that register made concrete: it triages **every** note across all 86 per-unit reports into one of four dispositions (fixed-now / feasible / not-possible / by-design) via the [[validation-findings-disposition]] process. Its most striking output is **green-washing detection** — a large share of "green" tests were found to assert a defective spec (tautologies, code-echo oracles), and were re-anchored to independently sourced literals and mutation-checked. Note the register is a 2026-06-12 snapshot **superseded by a full re-validation reset on 2026-06-24**, so read it as historical reasoning, not live status — the live per-unit status board is the [[validation-ledger]] (three phases: 86 implemented + 24 new campaign units, 148 Phase-2, 12 enhanced), which every governance page treats as ground truth for *where things stand*.

The campaign is tracked unit-by-unit in the [[test-unit-registry]] (364 units, 255 complete at the ingested revision), where each unit must clear the [[definition-of-done]]. Note the "258 algorithm units" figure from the README and the "364 test units" figure from [[algorithms-checklist-v2]] count different things — algorithm implementations versus tracked test units (255 done + 109 proposed) — and are not in conflict.
