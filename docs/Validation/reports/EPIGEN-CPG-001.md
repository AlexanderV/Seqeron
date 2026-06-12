# Validation Report: EPIGEN-CPG-001 — CpG Site / CpG Island Detection

- **Validated:** 2026-06-12   **Area:** Epigenetics
- **Canonical method(s):** `EpigeneticsAnalyzer.FindCpGSites(string)`, `EpigeneticsAnalyzer.CalculateCpGObservedExpected(string)`, `EpigeneticsAnalyzer.FindCpGIslands(string, minLength, minGc, minCpGRatio)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "CpG site"** (fetched 2026-06-12) confirmed:
  - A CpG site is `5'—C—phosphate—G—3'` — a **dinucleotide on a single strand**, explicitly *not* a C-G base pair.
  - Gardiner-Garden & Frommer (1987) CpG island criteria: length **≥ 200 bp**, GC content **> 50%**, observed/expected CpG ratio **> 0.6 (60%)**.
  - O/E formula: Observed CpGs / Expected, where **Expected = (number of C × number of G) / sequence length**. (Saxonov alternative: `((C+G)/2)² / L`.)
  - Takai & Jones (2002) stricter criteria: **> 500 bp, GC > 55%, O/E ≥ 0.65**.
- The TestSpec (`tests/TestSpecs/EPIGEN-CPG-001.md`) and Evidence doc (`docs/Evidence/EPIGEN-CPG-001-Evidence.md`) cite these same sources and numbers.

### Formula check
Canonical Gardiner-Garden & Frommer O/E:

```
O/E = CpG_count / ( (C_count × G_count) / N )      where N = window length
```

Equivalently `O/E = CpG_count × N / (C_count × G_count)` — the "×N factor" is present (it is folded into the denominator `(C×G)/N`). The spec uses Gardiner-Garden, **not** the Saxonov `((C+G)/2)²/L` variant. Confirmed correct.

### Edge-case semantics
- No CpG → O/E = 0 (numerator 0).
- C=0 or G=0 → expected = 0 → guarded return 0 (avoids div-by-zero).
- Length < 2 / null / empty → 0 (no dinucleotide possible).
- GpC ≠ CpG (only 5'→3' C-then-G counted).
- Adjacent CpGs in "CGCG" = 2 sites at 0 and 2 (distinct dinucleotide windows).
All have defined, sourced expected behaviour.

### Independent cross-check (hand computation)
| Seq | N | C | G | CpG | Expected = C·G/N | O/E | GC% |
|-----|---|---|---|-----|------------------|-----|-----|
| CGCGCGCGCGCGCGCGCGCG | 20 | 10 | 10 | 10 | 100/20 = 5.0 | **2.0** | 100% |
| ACGTCGACG | 9 | 3 | 3 | 3 | 9/9 = 1.0 | **3.0** | 67% |
| ACGT | 4 | 1 | 1 | 1 | 1/4 = 0.25 | **4.0** | 50% |
| AATT… (20) | 20 | 0 | 0 | 0 | 0 (guarded) | **0.0** | 0% |
| 400 bp "CGCG" | 400 | 200 | 200 | 200 | 40000/400 = 100 | **2.0** | 100% → island |

Worked island example (400 bp CGCG): length 400 ≥ 200, GC 100% > 50%, O/E 2.0 > 0.6 → **qualifies as CpG island** (matches M15). All spec O/E values and island calls reproduced independently.

### Findings / divergences
None. Description is biologically and mathematically correct.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs`
- `FindCpGSites` lines 114–128
- `CalculateCpGObservedExpected` lines 192–211
- `FindCpGIslands` lines 216–274 (+ `CalculateGcContent` 276–277)

### Formula realised correctly?
- `FindCpGSites`: uppercases input, scans `i = 0 .. Length-2`, yields `i` where `seq[i]=='C' && seq[i+1]=='G'` — exact 5'→3' C-then-G dinucleotide, GpC not matched. INV-1, INV-5 hold.
- `CalculateCpGObservedExpected` (line 209–210): `expected = (c * g) / (double)Length; return expected > 0 ? cpg / expected : 0;` — this is **exactly** Gardiner-Garden `CpG / ((C×G)/N)` with the `expected > 0` guard covering C=0 or G=0 (INV-2, INV-3). Length<2/null guarded at line 194. Note: `c * g` is `int` multiplication; for a 200–500 bp window max C·G ≈ (N/2)² well within `int` range — no overflow on stated ranges.
- `FindCpGIslands`: sliding window, step 1, window length `minLength` (the `Math.Min` at line 233 never shortens within the loop bound `i ≤ Length-minLength`, so every scanned window = `minLength`). Windows passing `gc >= minGc && cpgRatio >= minCpGRatio` are merged; the merged island is **re-validated** against GC% and O/E (lines 252–254, 267–269) and only emitted if `End-Start ≥ minLength`. Uses `>=` operators, matching Takai & Jones' explicit "≥ 200 bp, ObsCpG/ExpCpG ≥ 0.6, %GC ≥ 50%" phrasing. INV-4 holds.

### Cross-verification table recomputed vs code
All five datasets above were executed via the test fixture; M9→2.0, M11→3.0, M12→4.0, M10→0.0, C1 (CCCC)→0.0, M15→(0,400,Gc=1.0,O/E=2.0). Every value matches the external reference.

### Variant/delegate consistency
`FindCpGIslands` delegates O/E and GC computation to the same `CalculateCpGObservedExpected` / `CalculateGcContent` used standalone — consistent. `FindMethylationSites` (CpG context detection) uses the same C-then-G rule.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CpGDetection_Tests.cs` (24 tests). Assertions check exact sourced values (`Is.EqualTo(2.0).Within(1e-10)`, exact positions, exact Start/End/GcContent/CpGRatio), not tautologies; deterministic; cover all Stage-A edge cases (null, empty, single char, no C/G, GpC, adjacent, lowercase, thresholds, too-short).

### Findings / defects
None. Note (non-defect): spec §2 / checklist name the method `CalculateCpGObservedExpected` vs an older `CalculateObservedExpectedCpG` alias in the checklist index; the implemented and tested name is `CalculateCpGObservedExpected` — consistent across code, tests, and TestSpec §2.

## Verdict & follow-ups
- Stage A: **PASS** — definition, Gardiner-Garden O/E formula (denominator `(C×G)/N`), thresholds, and Takai-Jones alternative all confirmed against authoritative sources.
- Stage B: **PASS** — code realises the validated formula and thresholds exactly; div-by-zero guarded; all worked examples reproduced.
- **State: CLEAN.** No defect found.
- Build: succeeded. Tests: 29 CpG-filtered passed; full suite **4486 passed, 0 failed** (baseline preserved). No code changed.
