# Mutation Testing Analysis Report

**Date:** 2026-02-14  
**Tool:** Stryker.NET 4.12.0  
**Scope:** 6 projects (31 files with detected mutants), heavy algorithms excluded  
**Total mutants:** 3,804 generated | **1,327 survived** | **2,037 killed** | **440 NoCoverage**

---

## Summary

| Metric | Value |
|--------|-------|
| Overall mutation score | **60.6%** (killed / (killed + survived)) |
| Files analyzed | 31 (with active mutants) |
| CRITICAL files (< 50%) | 1 |
| HIGH priority files (50-60%) | 9 |
| MEDIUM priority (60-70%) | 7 |
| LOW priority (70%+) | 14 |

---

## Mutation Scores by File

| # | File | Score | Killed | Survived | NoCov | Total | Priority |
|---|------|-------|--------|----------|-------|-------|----------|
| 1 | CodonOptimizer.cs | 43.1% | 56 | 74 | 13 | 182 | CRITICAL |
| 2 | SpliceSitePredictor.cs | 51.2% | 128 | 122 | 0 | 317 | HIGH |
| 3 | PrimerDesigner.cs | 53.7% | 102 | 88 | 9 | 227 | HIGH |
| 4 | BedParser.cs | 54.4% | 92 | 77 | 14 | 224 | HIGH |
| 5 | ProteinSequence.cs | 54.9% | 50 | 41 | 0 | 119 | HIGH |
| 6 | GenBankParser.cs | 55.4% | 56 | 45 | 40 | 204 | HIGH |
| 7 | EmblParser.cs | 56.2% | 59 | 46 | 19 | 164 | HIGH |
| 8 | CodonUsageAnalyzer.cs | 56.3% | 40 | 31 | 9 | 99 | HIGH |
| 9 | ProbeDesigner.cs | 58.6% | 89 | 63 | 9 | 245 | HIGH |
| 10 | GcSkewCalculator.cs | 59.7% | 37 | 25 | 1 | 78 | HIGH |
| 11 | SequenceIO.cs | 60.6% | 114 | 74 | 19 | 279 | MEDIUM |
| 12 | EpigeneticsAnalyzer.cs | 61.9% | 122 | 75 | 19 | 268 | MEDIUM |
| 13 | FastqParser.cs | 62.9% | 73 | 43 | 1 | 151 | MEDIUM |
| 14 | DisorderPredictor.cs | 63.4% | 64 | 37 | 30 | 164 | MEDIUM |
| 15 | VcfParser.cs | 64.4% | 96 | 53 | 33 | 240 | MEDIUM |
| 16 | RepeatFinder.cs | 69.7% | 101 | 44 | 2 | 184 | MEDIUM |
| 17 | GffParser.cs | 69.9% | 51 | 22 | 8 | 113 | MEDIUM |
| 18 | SequenceComplexity.cs | 70.0% | 77 | 33 | 3 | 140 | LOW |
| 19 | MiRnaAnalyzer.cs | 70.6% | 151 | 63 | 1 | 254 | LOW |
| 20 | PopulationGeneticsAnalyzer.cs | 71.1% | 297 | 121 | 11 | 488 | LOW |
| 21 | GenomeAnnotator.cs | 71.5% | 98 | 39 | 12 | 209 | LOW |
| 22 | CrisprDesigner.cs | 72.3% | 73 | 28 | 7 | 152 | LOW |
| 23 | Translator.cs | 77.1% | 27 | 8 | 8 | 58 | LOW |
| 24 | FastaParser.cs | 77.4% | 24 | 7 | 3 | 46 | LOW |
| 25 | MotifFinder.cs | 78.1% | 82 | 23 | 26 | 169 | LOW |
| 26 | KmerAnalyzer.cs | 79.2% | 61 | 16 | 2 | 105 | LOW |
| 27 | RestrictionAnalyzer.cs | 81.7% | 49 | 11 | 0 | 82 | LOW |
| 28 | GenomicAnalyzer.cs | 82.6% | 57 | 12 | 16 | 116 | LOW |
| 29 | RnaSequence.cs | 91.9% | 34 | 3 | 0 | 47 | LOW |
| 30 | GeneticCode.cs | 92.3% | 12 | 1 | 0 | 18 | LOW |
| 31 | DnaSequence.cs | 92.6% | 25 | 2 | 0 | 38 | LOW |

**100% score (fully killed):** IupacHelper.cs (22/22), SequenceExtensions.cs (fixed during this session)

---

## Survived Mutants by Type

| Mutator Type | Count | % | Root Cause |
|-------------|-------|---|------------|
| Equality mutation (`<` ↔ `<=`, `>` ↔ `>=`, `==` ↔ `!=`) | 714 | 53.8% | Missing boundary-value tests |
| Arithmetic mutation (`+` ↔ `-`, `*` ↔ `/`) | 330 | 24.9% | Imprecise assertions or untested formulas |
| Logical mutation (`&&` ↔ `\|\|`, `!` negation) | 161 | 12.1% | Compound conditions not branch-tested |
| Block removal mutation | 74 | 5.6% | Dead or untested code paths |
| Null coalescing mutation (all variants) | 40 | 3.0% | No null-argument tests |
| Bitwise mutation (`&` ↔ `\|`) | 8 | 0.6% | Missing bit-pattern tests |

---

## Weakness Pattern Analysis

### Pattern 1: Boundary Conditions (54% of all survivors)

The **dominant weakness** across the entire codebase. Stryker swaps `<` to `<=` (or vice versa) and tests pass because they don't exercise exact boundary values.

**Examples:**
- `SpliceSitePredictor.cs:352` — `pos >= 0 && pos < sequence.Length` → mutated to `pos >= 0 || pos < sequence.Length` (survives: no test with pos=-1)
- `PrimerDesigner.cs:58` — `len <= param.MaxLength` → `len < param.MaxLength` (survives: no test with len == MaxLength)
- `CodonOptimizer.cs:347` — `currentGc >= minGc && currentGc <= maxGc` → mutated `>=` to `>` (survives: no test with gc == minGc exactly)

**Fix strategy:** For every boundary condition `if (x < N)`, add tests where x = N-1, N, and N+1.

---

### Pattern 2: Arithmetic Formula Precision (25%)

Mathematical formulas in scientific algorithms have mutations like `+` → `-`, `*` → `/`, constant changes. Tests use approximate comparison (`Assert.That(result, Is.EqualTo(x).Within(tolerance))`) or inputs where formula changes don't alter the outcome.

**Examples:**
- `PopulationGeneticsAnalyzer.cs:674` — `sumProduct / n + 2 * p1 * p2` → various arithmetic swaps survive
- `EpigeneticsAnalyzer.cs:501` — `(d - mean) * (d - mean) * (diffs.Count - 1)` → `*` ↔ `/`
- `MiRnaAnalyzer.cs:525` — `-stemLength * 1.5 - loopSize * 0.5` → sign and coefficient changes

**Fix strategy:** Use analytically computed expected values. Test with inputs where every arithmetic operator change yields a detectably different result (e.g., asymmetric inputs that break commutativity).

---

### Pattern 3: Boolean Logic Compound Conditions (12%)

`&&` ↔ `||` swaps survive in multi-clause conditions because tests satisfy ALL clauses or NONE — never hitting the case where only one clause is false.

**Examples:**
- `PrimerDesigner.cs:145` — `tm < param.MinTm && tm > param.MaxTm` → swap to `||` (should be `||` already?)
- `VcfParser.cs:775` — `sampleNames.Length > 0 && record.Format != null || record.Samples != null`
- `BedParser.cs:387` — `b.ChromEnd <= start && b.ChromStart >= end` → swap survives

**Fix strategy:** Create truth-table tests: for each clause, provide input where only THAT clause is false while others are true.

---

### Pattern 4: Block Removal / Dead Code (6%)

Entire `if`-blocks or method branches can be removed by Stryker and all tests still pass. This means the code path is either:
- Not tested at all
- Has no observable effect in tests (side-effect-only code)

**Hot spots:**
- `SequenceIO.cs` (7 block removals) — feature parsing branches in GenBank/EMBL parsers
- `MiRnaAnalyzer.cs` (6 block removals) — secondary analysis branches
- `ProbeDesigner.cs` (8 block removals) — validation sub-checks

**Fix strategy:** Add tests that specifically exercise the removed branch and assert its effect.

---

### Pattern 5: Null Handling (3%)

`??` operator mutations (remove left side, swap sides) survive because tests never provide `null` values for nullable parameters.

**Hot spots:** BedParser.cs (11 null coalescing), SequenceIO.cs (7), VcfParser.cs (4)

**Fix strategy:** Add `null`-argument tests for nullable parameters in parsers.

---

## Top-10 Worst Files — Hot Line Details

### 1. CodonOptimizer.cs — 43.1% (CRITICAL)

| Lines | Mutations | Pattern |
|-------|-----------|---------|
| L557-558 | 15x | RNA base-pair matching: `(b1 == 'A' && b2 == 'U') \|\| ...` — complex boolean with 6 clauses |
| L347 | 5x | GC content range check: `currentGc >= minGc \|\| currentGc <= maxGc` |
| L353, L497 | 10x | Loop boundary: `i > codons.Count` / `i > codons.Count - windowSize / 3` |
| L491 | 3x | Input validation: `IsNullOrEmpty && Length < windowSize` |
| L551-552 | 6x | Secondary structure scoring formula |

### 2. SpliceSitePredictor.cs — 51.2% (HIGH)

| Lines | Mutations | Pattern |
|-------|-----------|---------|
| L352, L388 | 9x | Position bounds checking: `pos >= 0 && pos < length` |
| L688 | 5x | Filtering: `Length < 500 \|\| Score < 0.8` |
| L434 | 4x | Score arithmetic: `donor + acceptor - branch` |
| L662, L670 | 6x | Position scaling: `Position * 50` |

### 3. PrimerDesigner.cs — 53.7% (HIGH)

| Lines | Mutations | Pattern |
|-------|-----------|---------|
| L370 | 8x | Hairpin loop detection: stem/loop length boundary |
| L58, L73 | 12x | Forward/reverse primer length constraints |
| L145 | 5x | Tm range check (suspected inverted logic: `&&` should be `\|\|`) |
| L299 | 4x | Input validation compound condition |

### 4-10: Similar patterns of boundary conditions and formula mutations.

---

## ~~Suspected Production Bugs~~ → FALSE POSITIVES (Test Weaknesses)

> **CORRECTION (verified 2025-02-14):** All 6 items below were **false positives** caused by
> misinterpreting Stryker's `Replacement` field. The `Replacement` shows the **mutated** code
> (what Stryker changed it TO), not the original. The production code already correctly uses `||`.
> These are **weak test coverage** issues, not production bugs.
>
> **Resolution:** Added 17 mutation-killing boundary tests across 4 test files.

The following `||→&&` logical mutations survived because tests did not exercise
individual clauses in isolation (only tested cases where both or neither clause was true):

1. **PrimerDesigner.cs:145** — Original: `tm < param.MinTm || tm > param.MaxTm` ✅ CORRECT  
   Mutation: `||→&&`. Fixed by: `EvaluatePrimer_TmOnlyBelowMin_FlagsTmIssue` + 3 boundary tests.

2. **PrimerDesigner.cs:299** — Original: `string.IsNullOrEmpty(sequence) || sequence.Length < ...` ✅ CORRECT  
   Mutation: `||→&&`. Fixed by: `HasHairpinPotential_NullSequence_ReturnsFalse` + 2 boundary tests.

3. **ProbeDesigner.cs:442** — Original: `gc < 0.40 || gc > 0.60` ✅ CORRECT  
   Mutation: `||→&&`. Fixed by: `DesignMolecularBeacon_AtRichTarget_ScorePenalizedForGcAndTm`.

4. **ProbeDesigner.cs:443** — Original: `tm < 55 || tm > 65` ✅ CORRECT  
   Mutation: `||→&&`. Fixed by: `DesignMolecularBeacon_GcRichTarget_ScorePenalizedForGcAndTm`.

5. **BedParser.cs:387** — Original: `b.ChromEnd <= start || b.ChromStart >= end` ✅ CORRECT  
   Mutation: `||→&&`. Fixed by: `Subtract_FragmentNotOverlappingWithSecondB_PreservedCorrectly` + 1 symmetric test.

6. **VcfParser.cs:568** — Original: `record.Samples == null || sampleIndex >= record.Samples.Count` ✅ CORRECT  
   Mutation: `||→&&`. Fixed by: `GetAlleleDepth_NullSamples_ReturnsNull` + 4 boundary tests.

### Algorithm Sources Consulted
- **PCR primer Tm ranges:** Wikipedia (Primer biology), ThermoFisher (Primer Analyzer) → MinTm=55, MaxTm=65°C standard
- **Probe GC/Tm criteria:** Standard molecular beacon design → GC 40-60%, Tm 55-65°C
- **BED interval arithmetic:** UCSC Genome Browser BED spec → 0-based half-open [chromStart, chromEnd)
- **VCF sample indexing:** VCF 4.3 specification → AD field, sample bounds checking


---

## Excluded Heavy Algorithms (not tested)

| File | LOC | Reason |
|------|-----|--------|
| SequenceAligner.cs | 976 | O(n×m) NW/SW, O(n²×k) MSA |
| ApproximateMatcher.cs | 415 | O(n×m) edit distance |
| RnaSecondaryStructure.cs | 845 | O(n³) Nussinov/MFE |
| PhylogeneticAnalyzer.cs | 768 | O(n²) distance matrix, UPGMA/NJ |
| ProteinMotifFinder.cs | 1,022 | Heavy domain prediction scans |
| ChromosomeAnalyzer.cs | 889 | Synteny analysis |
| MetagenomicsAnalyzer.cs | 701 | Genome binning + permutations |

These require either targeted mutation (specific methods only) or increased timeout settings.

---

## Recommendations

### ✅ Completed
1. **~~Исправить 6 вероятных багов~~** → All 6 were false positives. Added 17 mutation-killing boundary tests.

### Immediate Actions (убить максимум мутантов с минимумом тестов)
2. **CodonOptimizer.cs** — добавить boundary-тесты для GC-range и base-pair matching (15 мутантов на L557-558)
3. **Парсеры (BedParser, VcfParser, GenBankParser, EmblParser)** — добавить null-argument и boundary тесты

### Medium-term
4. Для всех файлов с score < 60% — систематично добавить boundary-value tests на каждое условие
5. Добавить truth-table тесты для compound boolean conditions (каждый clause отдельно)
6. Для формул в PopulationGeneticsAnalyzer и EpigeneticsAnalyzer — тесты с аналитически рассчитанными значениями

### Target
- **SHORT-TERM:** довести все CRITICAL/HIGH файлы до >= 70%
- **LONG-TERM:** overall mutation score >= 80%

---

## Batch Run Summary

| Project | Score | Duration | Files |
|---------|-------|----------|-------|
| Core | 73.3% | 13 min | 6 |
| Analysis | 64.2% | 7 min | 7 |
| MolTools | 57.0% | 53 min | 6 |
| Annotation | 60.2% | 22 min | 4 |
| Population | 69.2% | 5 min | 1 |
| IO | 52.9% | 6 min | 8 |
| **Total** | **60.6%** | **1h 47m** | **32** |
