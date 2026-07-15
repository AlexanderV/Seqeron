---
type: source
title: "Algorithms Checklist v2.0 — test-unit validation registry"
tags: [validation, testing]
doc_path: ALGORITHMS_CHECKLIST_V2.md
source_commit: 6a14170477c9472c0be07e3b7c7f7123e31eddcf
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Algorithms Checklist v2.0 — test-unit validation registry

The master QA ledger for Seqeron.Genomics (`ALGORITHMS_CHECKLIST_V2.md`, version 2.5, dated 2026-02-12): it enumerates every algorithm as a **test unit** with a stable ID, tracks its completion status, and records the acceptance criteria, evidence, and test artifacts behind it. It is the concrete backing for the [[validation-and-testing]] strategy.

## What it tracks

- **364 test units** total — **255 completed** (`☑`), 0 in progress/blocked, and **109 proposed / not started** (the assembly, phylogenetics, oncology, and cross-domain roadmap, prioritized in `ALGORITHMS_ROADMAP.md`).
- A **Processing Registry** table: per unit — area, method count, **evidence** (Wikipedia, primary literature by author-year, and reference tools such as Biopython, Rosalind, EMBOSS, Primer3, Kraken, TargetScan, REBASE), a `TestSpecs/<ID>.md`, and the C# test file(s).
- Per-area **Test Unit specs** giving each unit's **canonical method, complexity, invariant, and edge cases** — see [[test-unit-registry]] for the ID scheme and structure.
- The **[[definition-of-done]]** — the acceptance bar every unit must clear.
- Appendices: A (method index), B (complexity notes), C (canonical-implementation delegate trees), D (class coverage — **44/57 classes, 77%**).

## Honest caveats it records

- **Unverified complexity claims** (Appendix B): several units carry a claimed complexity marked `⚠️` that does not hold as stated — e.g. `REP-STR-001` is actually O(n × U × R), `META-BIN-001` is O(n × k × i), and `CRISPR-OFF-001` "may be exponential with high mismatches."
- **Campaign-added units pending re-validation**: 21 net-new algorithms from the limitation-elimination campaign are marked `☐` — the code ships and is covered by a test fixture, but the unit has **not** been re-validated under the project's two-stage (Stage A/B) protocol. These reinforce [[research-grade-limitations]].

## Where this fits

- [[validation-and-testing]] — this registry is the concrete instance of that strategy.
- [[test-unit-registry]] — the Test Unit ID scheme and per-unit structure defined here.
- [[definition-of-done]] — the acceptance criteria defined here.
- [[research-grade-limitations]] — the unverified-complexity and pending-re-validation caveats.
