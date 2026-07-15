---
type: concept
title: "Validation findings disposition (A/B/C/D triage)"
tags: [validation, testing, governance]
sources:
  - docs/Validation/FINDINGS_REGISTER.md
source_commit: 9710ae416a27c89fd4c0699c0a7631c0df408224
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: findings-register
      evidence: "The register collects every note/limitation/follow-up across all 86 per-unit validation reports in docs/Validation/reports/ and dispositions each — it is the finding-triage side of the validation campaign."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:build-quality-gate
      source: findings-register
      evidence: "Green-washing detection (re-derive the oracle from an external source, mutation-check the test) is the same 'tests must encode the algorithm/business spec, not a buggy impl' rule the Sonar gate applied to the S1244/S125 giants."
      confidence: medium
      status: current
---

# Validation findings disposition (A/B/C/D triage)

The governance process Seqeron uses to close out the notes its
[[validation-and-testing|validation campaign]] produces: **every** finding across all per-unit
reports is placed in **exactly one** of four disposition buckets, and each bucket carries a fixed
action policy. The canonical instance of this process is the [[findings-register]]
(`docs/Validation/FINDINGS_REGISTER.md`).

## The four buckets and their policies

| Bucket | Meaning | Action |
|--------|---------|--------|
| **FIXED-NOW** | Doc/comment/spec text or a small, safe code change. | Fix in the pass. |
| **FEASIBLE -> IMPLEMENT** | Needs real code but is safe with strict tests in its own context. | Implement separately, tests first. |
| **NOT-POSSIBLE (radical)** | Requires a redesign / public-API change / new model. | Document; defer to the ledger backlog. |
| **BY-DESIGN** | Correct, sourced intended behaviour. | No change; "fixing" it would be wrong. |

The discipline is that the categories are **mutually exclusive** — a finding is triaged once — and
that "radical" is a *documentation* verdict, not a permanent one: in the register most NOT-POSSIBLE
items were later reclassified and implemented as approved breaking or additive API changes once a
faithful, source-grounded design was found.

## Green-washing detection is the core mechanism

The most valuable output of the triage is not the code fixes but the discovery that many "green"
tests were **passing a defective spec** — tautologies, `Contains`/range assertions, or **code-echo**
oracles recomputed from the implementation's own constants. The standing remedy, applied uniformly:

1. Re-derive the expected value **independently** from an externally retrieved primary source
   (paper, reference implementation, or verbatim data file) — never from the code under test.
2. Lock the exact literal (tight tolerance).
3. **Mutation-check**: confirm a deliberately seeded bug now fails the strengthened test.

This is the "**tests must encode the algorithm/business spec, not a buggy impl**" rule — the same
principle the [[build-quality-gate]] used to reject blind linter autofixes on the S1244/S125 giants.
A finding that turns out to be a sourced non-defect is correctly closed as **BY-DESIGN**, not "fixed."

This remedy is not ad-hoc — it is baked into the [[validation-protocol]] itself: every unit is
re-validated in a **fresh session** separate from the implementer's, cross-checking each formula and
oracle against an **independent external primary source** before the code is trusted, which is exactly
what surfaces a code-echo oracle in the first place.

## Relationship to live status

A disposition register is a **snapshot** of one pass, not a live tracker. The findings register was
explicitly **superseded on 2026-06-24** by a full re-validation reset (all units back to pending in
the [[validation-ledger]], `docs/Validation/VALIDATION_LEDGER.md`); the per-unit [[test-unit-registry]] and its
[[definition-of-done]] carry current state. Treat a disposition table as historical evidence of *how*
findings were reasoned about, and the [[validation-ledger|ledger]] as ground truth for *where things stand*.
