---
type: concept
title: "Build-time quality gate (Sonar ratchet + warnings-as-errors)"
tags: [quality-gate, build]
sources:
  - docs/sonar-gate-plan.md
source_commit: b7f71a58c61a48497ab33b3a504187fa8562a392
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: sonar-gate-plan
      evidence: "Definition of done per rule requires `dotnet test` is green in addition to the blocking `dotnet build` gate — the static-analysis gate rides on top of the correctness test suite."
      confidence: high
      status: current
---

# Build-time quality gate (Sonar ratchet + warnings-as-errors)

Seqeron enforces a solution-wide **static-analysis quality gate at build time** — distinct
from runtime [[scientific-rigor]] and from the correctness-focused
[[validation-and-testing]] suite. `SonarAnalyzer.CSharp` is wired in via
`Directory.Build.props`, and the solution builds under `TreatWarningsAsErrors`: once a rule
is enforced, any new finding **fails the build**.

## The ratchet

Every firing rule begins **report-only**, parked in `WarningsNotAsErrors`. The gate is
tightened by *ratcheting* one rule per session from advisory to **blocking**, worked
most-important first (real bugs → test integrity → deliberate-silences → safe cleanups →
high-effort readability). A rule is resolved one of two ways:

- **Fix** — resolve every finding in code, then remove the rule's code from
  `WarningsNotAsErrors` so it trips the solution-wide gate.
- **Silence with justification** — set `dotnet_diagnostic.<code>.severity = none` in
  `.editorconfig` with a one-line rationale, and drop it from the list. Silences are reserved
  for rules that fight this ASCII, culture-free, performance-oriented bioinformatics library
  (non-crypto `Random`, `unsafe` for the memory-mapped suffix tree, `GC.Collect` in benchmarks,
  perf-critical accumulation loops, exposed numeric buffers).

The [[sonar-gate-plan]] tracker records the full run: **66 of 66 rules resolved**, zero
remaining Sonar warnings, green under `TreatWarningsAsErrors`.

## Review, not blind fixing

The gate's sharp edge is that a linter autofix can be *worse* than the finding. The two
highest-volume rules — exact float `==` (S1244, 234 sites) and "commented-out code" (S125,
834 blocks) — were resolved by **per-site review**, because applying Sonar's fix (add a
tolerance / delete the comment) would have changed numeric results or destroyed intentional
scientific-provenance documentation. This is the "**tests must encode the algorithm/business
spec, not a buggy impl**" standing rule applied to static analysis: verify against the
documented intent before "fixing." Silencing a reviewed non-defect is the correct outcome.

Behaviour changes made to satisfy a rule are treated as deliberate and test-verified — e.g.
splitting iterator argument-validation to **fail fast** (S4456) was run against the full suite
to confirm no test depended on the old lazy-throw behaviour.
