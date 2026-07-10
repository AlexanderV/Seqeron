---
type: gotcha
title: "Research-grade, not for clinical use"
tags: [limitations, status]
sources:
  - README.md
  - ALGORITHMS_CHECKLIST_V2.md
  - docs/Validation/LIMITATIONS.md
source_commit: 45545719fbdd7689c20bb680104862f6098adf32
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: constrains
      object: concept:scientific-rigor
      source: readme
      evidence: "Beta — research-grade software, not for clinical or diagnostic use ... Many algorithms are faithful but simplified or subset realisations of fuller published methods"
      confidence: high
      status: current
---

# Research-grade, not for clinical use

Seqeron is **beta** and explicitly research-grade: feature-complete with a stabilizing public API, but **not for clinical or diagnostic use**. This bounds how far the project's [[scientific-rigor]] and [[validation-and-testing]] claims can be relied upon.

## The sharp edges

- **No external validation** — no third-party/independent audit, peer review, or regulatory clearance. All validation to date is internal and self-conducted.
- **Simplified / subset implementations** — many algorithms are faithful but **simplified or subset** realisations of fuller published methods; their honest scope is documented in `docs/Validation/LIMITATIONS.md` and guarded at runtime by `LimitationPolicy` — see [[operating-envelope-and-limitation-policy]] for the per-unit envelope and its three enforcement modes.
- **Some units ship ahead of validation** — per [[algorithms-checklist-v2]], 21 campaign-added algorithms are marked `☐` pending independent Stage A/B re-validation (code present and fixture-covered, DoD not yet re-met), and 109 proposed units are not started.
- **Not all complexity claims are verified** — Appendix B of the checklist flags several unit complexities as `⚠️` (e.g. `CRISPR-OFF-001` "may be exponential with high mismatches"; `REP-STR-001` is really O(n × U × R)).
- **Unstable public API** — APIs may still change between releases on the road to 1.0.
- **No warranty** — the authors disclaim warranties of correctness or fitness; use at your own risk (see `LICENSE`).

## What this requires of a user

Before relying on any output with real data or in production, independently verify it against established tools for the specific use case, and never use it for clinical or diagnostic decision-making without separate qualification and validation. The full status discussion lives in the "Project status & validation" section of `README.md`.
