---
type: source
title: "Checklist 04: Mutation Testing (Stryker)"
tags: [validation, testing, methodology]
doc_path: docs/checklists/04_MUTATION_TESTING.md
sources:
  - docs/checklists/04_MUTATION_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 04: Mutation Testing (Stryker)

The **P1** per-unit checklist for mutation testing — a 258-row table tracking each unit to a
mutation score ≥ 80% target. Synthesized in the concept [[mutation-testing]] (which also covers
the earlier [[mutation-testing-analysis]] 2026-02-14 baseline); part of the
[[validation-and-testing]] program tracked in the [[test-unit-registry]].

## What this file records

- **Framework:** Stryker.NET; process = run `dotnet stryker` per source file, analyse survivors,
  write targeted killer-tests for survivors below 80%. Target **≥ 80% per file**, ~25 source
  files to mutate.
- **Per-unit table:** all **258 rows ☑** (mutation score ≥ 80%).
- **Ceiling fully broken (2026-06-30):** the last-mile residual on three hard files (deep DP /
  iterative-solver pipelines) was the **equivalent-mutant-with-respect-to-the-public-API** class
  — mutations that don't change *observable* output because the robust solver / re-convergent DP
  hides the internal arithmetic. Faithful resolution (explicitly **not** brittle characterization
  / green-washing): make the internal computation observable (`internal` + `InternalsVisibleTo`)
  and assert it against an **independent first-principles oracle**.

## Deviations and contradictions

**Temporal progression, not contradiction, vs [[mutation-testing-analysis]]:** the standalone
Evidence report scored **60.6%** overall on 2026-02-14 over the tractable surface; this checklist
records the later **all-files-≥-80%, 258/258** end-state reached 2026-06-30. Same methodology,
two points in time.
