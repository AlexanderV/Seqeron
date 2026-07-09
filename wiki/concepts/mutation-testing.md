---
type: concept
title: "Mutation testing (Stryker.NET) — measuring test-suite strength"
tags: [testing, validation, methodology]
sources:
  - docs/Evidence/MUTATION-TESTING-ANALYSIS.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: mutation-testing-analysis
      evidence: "Mutation (Stryker.NET) — weak tests, by seeding deliberate bugs — is one of the ten complementary methodologies; this run seeds 3,804 mutants to score the suite."
      confidence: high
      status: current
---

# Mutation testing (Stryker.NET) — measuring test-suite strength

Mutation testing measures **how good the tests are**, not whether the code is correct. A tool ([Stryker.NET](https://stryker-mutator.io/), 4.12.0 here) seeds thousands of small deliberate faults ("mutants") into the source — swapping `<` for `<=`, `+` for `-`, `&&` for `||`, removing a block, dropping a `??` — recompiles, and reruns the suite against each. A mutant that makes at least one test fail is **killed**; one that leaves the suite green **survives**, exposing a gap. The headline number is the **mutation score** = killed / (killed + survived). This is the [[validation-and-testing]] methodology that the [[advanced-testing-checklist]] rated ★★★★☆ / P1 — the deliberate-bug-seeding complement to [[property-based-testing]] and [[metamorphic-testing]]. The underlying baseline report is [[mutation-testing-analysis]] (60.6% on 2026-02-14); the per-unit tracker [[mutation-testing-checklist]] records the later end-state where **every source file reached ≥ 80%** (ceiling broken 2026-06-30, resolving the equivalent-mutant-with-respect-to-the-public-API residual via `InternalsVisibleTo` + an independent first-principles oracle rather than green-washed characterization tests).

## The 2026-02-14 baseline run

Scope: 6 projects, 31 files with active mutants, heavy algorithms excluded. 3,804 mutants generated → 2,037 killed, 1,327 survived, 440 no-coverage → **overall score 60.6%**. Per-project scores ranged from IO 52.9% and MolTools 57.0% up to Core 73.3%. File scores spanned `CodonOptimizer.cs` 43.1% (the single CRITICAL file, < 50%) to value objects `DnaSequence`/`RnaSequence`/`GeneticCode` at 92%+ and `IupacHelper.cs` at 100%. Targets set: bring all CRITICAL/HIGH files to ≥ 70% short-term, overall to ≥ 80% long-term.

## Survivor taxonomy — where suites are weak

The surviving mutants cluster into five reusable weakness patterns (percentages of all survivors):

- **Boundary conditions (54%)** — the dominant weakness. `<`↔`<=` / `>`↔`>=` swaps survive because tests never exercise the exact boundary value. Fix: for every `if (x < N)`, test x = N−1, N, N+1.
- **Arithmetic formula precision (25%)** — `+`↔`-`, `*`↔`/`, constant changes in scientific formulas survive because tests use loose tolerance or symmetric inputs. Fix: analytically computed expected values; asymmetric inputs that break commutativity.
- **Boolean compound conditions (12%)** — `&&`↔`||` swaps survive because tests satisfy all clauses or none, never exactly one clause false. Fix: truth-table tests, one input per clause-false.
- **Block removal / dead code (6%)** — whole `if`-blocks deletable with the suite still green ⇒ path untested or side-effect-only. Fix: exercise the branch and assert its effect.
- **Null handling (3%)** — `??` mutations survive because nullable parameters are never passed `null`. Fix: null-argument tests (parsers are the hot spot).

## The false-positive lesson

Six survivors were first logged as *suspected production bugs* (`||` looked like it should be `&&`), then retracted: Stryker's `Replacement` field shows the **mutated** code, not the original. The production code was already correct; these were weak-coverage issues, not defects. Resolution was 17 mutation-killing boundary tests, not code fixes — a concrete instance of the standing rule that [[build-quality-gate|tests must encode the spec, not chase the tool's output]]. Fixed expected values were cross-checked against external algorithm sources (primer Tm 55–65 °C, molecular-beacon GC 40–60%, UCSC BED half-open intervals, VCF 4.3 AD indexing).

## Scope caveat

Seven heavy algorithms (`SequenceAligner`, `RnaSecondaryStructure` O(n³), `PhylogeneticAnalyzer`, `ProteinMotifFinder`, `ChromosomeAnalyzer`, `MetagenomicsAnalyzer`, `ApproximateMatcher`) were **excluded** — their O(n·m)/O(n³) cost makes full mutation runs infeasible without targeted per-method mutation or raised timeouts. So 60.6% is the score over the *tractable* surface, not the whole library — one more entry in the [[research-grade-limitations|research-grade caveat]] that the testing program is deep but has uncovered surface.
