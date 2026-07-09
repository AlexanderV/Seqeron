---
type: concept
title: "Differential testing — cross-checking two independent implementations"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/08_DIFFERENTIAL_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: differential-testing-checklist
      evidence: "Priority P2. Compares the outputs of two independent implementations of the same algorithm on identical inputs to surface subtle implementation bugs — reference implementations (Biopython, EMBOSS) and alternative algorithmic strategies serve as the second oracle."
      confidence: high
      status: current
---

# Differential testing

Differential testing runs **two independent implementations** of the same algorithm on identical
inputs and asserts they agree. Disagreement localises a subtle implementation bug that neither a
single-implementation property nor a hand-written assertion would catch — the second
implementation *is* the oracle. In bioinformatics the second implementation can be a reference
library (Biopython, EMBOSS), a brute-force implementation valid on small inputs, or an
alternative algorithmic strategy inside the project. This is a **P2** member of the
[[validation-and-testing]] program; the checklist record is [[differential-testing-checklist]].

## The strategy taxonomy

- **ALT — Alternative algorithm**: a different approach to the same problem cross-checks the
  production one.
- **BRUTE — Brute force**: an exhaustive reference implementation, correct-by-construction on
  small inputs, checks the optimised path. Example: **PROTMOTIF-HMM-001** cross-checks the Plan7
  local+glocal DP against an independent path-enumeration (`Plan7ProfileHmm_
  ForwardBackwardDifferential_Tests`).
- **REF — Reference comparison**: against a hand-computed result or an external reference library.
- **DUAL — Dual in-project implementations**: two methods of one class computing the same thing
  by different routes.

## Coverage

**107 / 255 complete, 151 not started** at the 2026-03-19 checklist — a partially-realised
methodology, with the starting point being `SafeVsUnsafeDifferentialTests` for SuffixTree plus
the Plan7 HMM differential. The checklist ranks the remaining work by feasibility: **~25**
high-value pairs where an ALT or BRUTE reference is practical, **~35** medium (REF comparison),
**~26** lower (needing a DUAL re-implementation). Differential testing is the most expensive of
the ten (it requires a *second* correct implementation), which is why it — like
[[snapshot-testing]] — is deliberately incomplete rather than force-fit, in contrast to the
258/258 invariant methodologies ([[property-based-testing]], [[metamorphic-testing]]).
