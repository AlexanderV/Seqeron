---
type: source
title: "Evidence: Mutation Testing Analysis (Stryker.NET 2026-02-14 baseline)"
tags: [validation, testing, methodology]
doc_path: docs/Evidence/MUTATION-TESTING-ANALYSIS.md
sources:
  - docs/Evidence/MUTATION-TESTING-ANALYSIS.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: Mutation Testing Analysis

The Stryker.NET mutation-testing report that scores the **strength of the test suite**
(not the correctness of the code) by seeding deliberate faults and checking which the
tests catch. The baseline run, survivor taxonomy, the false-positive lesson, and the
excluded-algorithm scope caveat are synthesized in the concept [[mutation-testing]]; this
is one instance of the [[algorithm-validation-evidence|evidence artifact]] pattern and
part of the broader [[validation-and-testing]] program. Unlike the per-algorithm Evidence
files it is a **methodology report**, not a [[test-unit-registry]] unit — it measures the
suite that validates those units.

## What this file records

- **Run metadata:** Stryker.NET 4.12.0, 2026-02-14; 6 projects, 31 files with active
  mutants, seven heavy algorithms excluded. 3,804 mutants → **2,037 killed, 1,327
  survived, 440 no-coverage** ⇒ **overall mutation score 60.6%** = killed / (killed +
  survived). Per-project: Core 73.3%, Analysis 64.2%, Annotation 60.2%, Population 69.2%,
  MolTools 57.0%, IO 52.9%.

- **File-level extremes:** `CodonOptimizer.cs` 43.1% (the single CRITICAL file, < 50%) up
  to value objects `DnaSequence`/`RnaSequence`/`GeneticCode` at 92%+ and `IupacHelper.cs`
  100%. Nine HIGH files (50–60%), seven MEDIUM (60–70%), fourteen LOW (70%+).

- **Survivor taxonomy (five reusable weakness patterns):** boundary conditions
  (`<`↔`<=`, `>`↔`>=`, `==`↔`!=`) **53.8%** — the dominant gap; arithmetic-formula
  precision (`+`↔`-`, `*`↔`/`) **24.9%**; boolean compound conditions (`&&`↔`||`, `!`)
  **12.1%**; block removal / dead code **5.6%**; null-coalescing (`??`) **3.0%**; bitwise
  0.6%. Each row carries a root cause and a fix strategy (e.g. for `if (x < N)` test
  x = N−1, N, N+1; truth-table tests for compound conditions; null-argument tests for
  parsers).

- **The false-positive lesson (recorded correction):** six survivors first logged as
  *suspected production bugs* (`||` looked like it should be `&&`) were **retracted** —
  Stryker's `Replacement` field shows the **mutated** code, not the original; the
  production code was already correct. Resolution was **17 mutation-killing boundary
  tests**, not code fixes. Fixed expected values were cross-checked against external
  algorithm sources (primer Tm 55–65 °C, molecular-beacon GC 40–60%, UCSC BED half-open
  intervals, VCF 4.3 AD indexing).

- **Scope caveat:** seven heavy algorithms (`SequenceAligner`, `RnaSecondaryStructure`
  O(n³), `PhylogeneticAnalyzer`, `ProteinMotifFinder`, `ChromosomeAnalyzer`,
  `MetagenomicsAnalyzer`, `ApproximateMatcher`) were **excluded** — their O(n·m)/O(n³)
  cost makes full mutation runs infeasible without targeted per-method mutation or raised
  timeouts. So 60.6% covers the *tractable* surface, not the whole library — one more
  entry in the [[research-grade-limitations|research-grade caveat]].

- **Targets set:** bring all CRITICAL/HIGH files to ≥ 70% short-term; overall to ≥ 80%
  long-term.

## Deviations and contradictions

**No source contradictions.** The one recorded reversal is internal and self-correcting:
the six "suspected production bugs" were reclassified as weak-coverage false positives
once the `Replacement`-field misreading was caught — a concrete instance of the standing
rule that tests must encode the spec, not chase the tool's output. The excluded heavy
algorithms are the only acknowledged coverage gap; the 60.6% figure is explicitly scoped
to the tractable surface, not presented as a whole-library score.
