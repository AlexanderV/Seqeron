---
type: source
title: "Sonar gate ratchet — work plan"
tags: [quality-gate, build, refactoring]
doc_path: docs/sonar-gate-plan.md
source_commit: b7f71a58c61a48497ab33b3a504187fa8562a392
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Sonar gate ratchet — work plan

The tracker for moving `SonarAnalyzer.CSharp` from advisory to enforced across the
solution: every firing rule starts *report-only* (listed in `WarningsNotAsErrors` in
`Directory.Build.props`), and this plan ratchets each to **blocking** one rule per
session, most-important first. It is the execution log for the [[build-quality-gate]]
concept, and the development-time companion to the correctness-focused
[[validation-and-testing]] strategy.

## Outcome (as of 2026-07-09)

**RATCHET COMPLETE — 66 of 66 rules resolved.** Every SonarAnalyzer rule that fired is
now either **blocking** (findings fixed in code) or **silenced with a documented
justification** in `.editorconfig` (`severity = none`). `Directory.Build.props` no longer
carries any report-only Sonar codes, so the whole solution builds green under
`TreatWarningsAsErrors` with **zero** remaining Sonar warnings. The full suite — **14
assemblies, 20,266 core tests** — is green (consistent with the 22k+ figure in
[[validation-and-testing]]).

## Definition of done, per rule

Each rule clears a five-step bar (see `docs/sonar-gate-plan.md` §"Definition of done"):
1. Every finding resolved — fixed in code, or the rule silenced with justification.
2. The rule's code removed from `WarningsNotAsErrors` in `Directory.Build.props` (so it
   now trips the solution-wide `TreatWarningsAsErrors` gate). For a *silence* item instead
   set `dotnet_diagnostic.<code>.severity = none` in `.editorconfig` and drop it from the list.
3. `dotnet build Seqeron.sln -c Debug` green (gate on).
4. `dotnet test` green.
5. Tick the checkbox.

The two resolution paths — **fix** vs **silence** — and the "one rule per session" cadence
are the heart of the ratchet; see [[build-quality-gate]].

## How the 66 rules were grouped

The plan bins rules by risk/value, worked highest-value first:

- **Group A — real bugs (fixed, gate):** genuine defects Sonar caught, e.g. **S2184**
  integer division truncating before a `double` assign in `MiRnaAnalyzer` (the one entry in
  the Log with a dedicated commit note), **S1871** identical branches in `SpliceSitePredictor`,
  **S1994/S127** loop stop-condition mismatches, **S2930** undisposed `IDisposable`, **S1854**
  dead stores, **S2696** instance method writing a static field.
- **Group B — test integrity (fixed, gate):** **S2699** tests with no assertions, **S2187**
  a test class with no tests.
- **Group C — deliberate/inapplicable (silenced):** rules that fight this ASCII, culture-free,
  performance-oriented bioinformatics library — non-crypto `Random` (**S2245**), `unsafe`
  blocks for the memory-mapped suffix tree (**S6640**), `GC.Collect` in benchmark paths
  (**S1215**), deprecated-code notices (**S1133**), interop/acronym type names (**S101**).
  **S2479** (literal control chars) was *fixed*, not silenced.
- **Group D — safe cleanups (mixed):** the bulk by volume — unused locals (**S1481**),
  `.Last()`/`[^1]` indexing (**S6608**), collapsible `if` (**S1066**), `params` array removal
  (**S3878**), `static`/`readonly`/`sealed` promotions, etc. Many were fixed only in production
  and *relaxed in tests* where the pattern is deliberate scaffolding (exercise-call locals,
  micro-perf indexing, `...Async().GetAwaiter().GetResult()`).
- **Group E — large / high-effort readability (last):** **S3358** nested ternary→`switch`,
  **S4456** iterator validation split, **S3267** loop→`Where` (silenced, hot-path), and the two
  reviewed giants **S1244** (234 float `==`) and **S125** (834 commented-code blocks).

## Notable behaviour changes and judgement calls

- **S4456 — intentional behaviour change.** All 27 iterator methods were split into an eager
  validation wrapper + a private `...Core` iterator, so **invalid arguments now fail fast
  instead of throwing on first enumeration.** Two methods validated *after* a `yield break`
  empty-guard (empty input masked a bad arg) were reordered to validate unconditionally. The
  full suite was run specifically to catch tests asserting the old lazy behaviour — none existed.
- **S1244 and S125 — resolved by review, not blind fixing.** Per the standing "tests must
  encode the algorithm/business spec, not a buggy impl" rule: all 13 production **S1244** float
  `==` sites are intentional numerical semantics (convergence `sum == previous`, zero-range /
  division guards, deterministic tie-breaks, degenerate-case handling) where a tolerance would be
  *wrong*; the 416 flagged **S125** lines are scientific-provenance documentation Sonar
  misclassifies as dead code (source-citation formulas, reference-impl snippets from TargetScan
  `.pl` / HMMER / ViennaRNA / Primer3, worked examples). Both were silenced as reviewed non-defects.

## Caveats and contradictions

- **Internal staleness inside the doc.** The file opens with the "**66 of 66 COMPLETE**"
  banner but *also* still contains an earlier "**Optimal finish sequence (the remaining 31
  rules)**" section and a "Remaining: 5 deferred Group D rules … and the 5 Group E giants"
  sentence. The Log table's final rows (2026-07-09: S3267/S1244/S125 and S3358; 2026-07-08:
  S4456) confirm completion, so the "remaining" planning prose is a **stale earlier draft**
  superseded by the completion banner — not a live backlog. Read the Log as ground truth.
- **Pre-existing flaky property tests, unrelated.** A few FsCheck property tests (e.g.
  `MultipleAlign_TotalScore_NonNegative_ForHomologs`, `NnTm_HigherSodium_NotLowerTm`) fail
  ~1 per full run on unlucky seeds and pass in isolation; the doc flags them as pre-existing and
  unrelated to the gate work. Reinforces the aspirational-coverage caveat in
  [[research-grade-limitations]].

## Where this fits

- [[build-quality-gate]] — the concept this document is the execution log for.
- [[validation-and-testing]] — the correctness-testing side of the same dev-time correctness story.
- [[research-grade-limitations]] — the flaky-property and reviewed-not-fixed notes touch its caveats.
