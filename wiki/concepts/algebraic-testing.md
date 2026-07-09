---
type: concept
title: "Algebraic testing — verifying explicit algebraic laws"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/06_ALGEBRAIC_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: algebraic-testing-checklist
      evidence: "Priority P1. Verifies algebraic laws (identity, commutativity, associativity, involution, idempotence, round-trip/isomorphism, distributivity) — 'більш формальний підхід, ніж загальні property-тести' (a more formal approach than general property tests)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:property-based-testing
      source: algebraic-testing-checklist
      evidence: "Algebraic laws exist implicitly in the FsCheck Property files (involution complement in SequenceProperties, round-trip in FastaRoundTripProperties) but algebraic testing systematises them as named laws; ~172 laws over 89 applicable units, the rest marked ✗ not-applicable."
      confidence: high
      status: current
---

# Algebraic testing

Algebraic testing verifies that operations obey **explicit algebraic laws** — identity,
commutativity, associativity, involution, idempotence, round-trip (isomorphism), distributivity,
triangle inequality. It is a **more formal, more targeted** relative of [[property-based-testing]]:
where a property asserts a loose invariant over generated inputs, an algebraic test names the
exact law an operation must satisfy. Many genomic operations have such laws by construction —
`complement(complement(x)) = x`, `parse(serialize(x)) = x`, `score(a,b) = score(b,a)`,
`d(x,x) = 0` with the triangle inequality on distances. This is a **P1** member of the
[[validation-and-testing]] program; the checklist record is [[algebraic-testing-checklist]].

## The law taxonomy

- **ID — Identity**: `f(x, e) = x` or a neutral element maps to zero.
- **COMM — Commutativity**: `f(a,b) = f(b,a)`.
- **ASSOC — Associativity**: `f(f(a,b),c) = f(a,f(b,c))`.
- **INV — Involution**: `f(f(x)) = x` (complement, reverse-complement).
- **IDEMP — Idempotence**: `f(f(x)) = f(x)` (normalisation, validation).
- **RT — Round-trip / Isomorphism**: `g(f(x)) = x` (parse ∘ serialize on every file format).
- **DIST — Distributivity / Conservation**: additivity and conservation laws (MW additive over
  residues minus water; fragment lengths sum to sequence length).
- **TRI — Triangle inequality**: on metric distances (edit, p-distance, Fst).

## Coverage

**89 units complete, ~172 laws** (≈2 per applicable algorithm); **169 units marked ✗ not
applicable** — the honest recognition that not every operation has an algebraic law, so the
methodology is deliberately partial rather than force-fit. The laws were already latent in the
FsCheck property files (involution in `SequenceProperties`, round-trip in
`FastaRoundTripProperties`); algebraic testing's contribution is to **systematise them as named
laws** rather than leave them implicit. It sits between [[property-based-testing]] (looser,
universal) and [[metamorphic-testing]] (relations across runs when no oracle exists), all three
tracked per unit in the [[test-unit-registry]].
