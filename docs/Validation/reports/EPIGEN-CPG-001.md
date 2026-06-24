# Validation Report: EPIGEN-CPG-001 — CpG Site / CpG Island Detection

- **Validated:** 2026-06-24   **Area:** Epigenetics
- **Canonical method(s):** `EpigeneticsAnalyzer.FindCpGSites(string)`, `EpigeneticsAnalyzer.CalculateCpGObservedExpected(string)`, `EpigeneticsAnalyzer.FindCpGIslands(string, minLength, minGc, minCpGRatio)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "CpG site"** (fetched 2026-06-24) confirmed:
  - A CpG site is `5'—C—phosphate—G—3'` — a dinucleotide on a single strand, explicitly *not* a C·G base pair.
  - Gardiner-Garden & Frommer (1987) CpG island criteria: length **≥ 200 bp**, GC content **> 50%**, observed/expected CpG ratio **> 0.6 (60%)**.
  - GGF O/E formula: Observed CpGs / Expected, where **Expected = (number of C × number of G) / length**. A documented *alternative* (Saxonov) is `((C+G)/2)² / length` — the spec uses the **GGF C×G** convention, not Saxonov.
  - Takai & Jones (2002) stricter criteria: **> 500 bp, GC > 55%, O/E ≥ 0.65**.
- **Web search** (peer-reviewed / tool refs) independently restated the same GGF formula in the algebraically equivalent form `O/E = [CpG / (C × G)] × N`, i.e. `R = (A×B)/(C×D)` with A=CpG, B=length, C=#C, D=#G.

### Formula check
Canonical Gardiner-Garden & Frommer O/E:

```
O/E = CpG_count / ( (C_count × G_count) / N )  =  CpG_count × N / (C_count × G_count)   (N = window length)
```

Confirmed against both Wikipedia and the independent search restatement. The denominator is **C × G** (GGF), not `((C+G)/2)²` (Saxonov). Confirmed correct.

### Edge-case semantics
- No CpG → O/E = 0 (numerator 0).
- C=0 or G=0 → expected = 0 → guarded return 0 (avoids div-by-zero).
- Length < 2 / null / empty → 0 (no dinucleotide possible).
- GpC ≠ CpG (only 5'→3' C-then-G counted).
- Adjacent CpGs in "CGCG" = 2 sites at 0 and 2.
All have defined, sourced expected behaviour.

### Independent cross-check (hand computation)
| Seq | N | C | G | CpG | Expected = C·G/N | O/E (GGF) | O/E (Saxonov) |
|-----|---|---|---|-----|------------------|-----------|---------------|
| CGCG×5 (20bp) | 20 | 10 | 10 | 10 | 100/20 = 5.0 | **2.0** | 2.0 |
| ACGTCGACG | 9 | 3 | 3 | 3 | 9/9 = 1.0 | **3.0** | 3.0 |
| ACGT | 4 | 1 | 1 | 1 | 1/4 = 0.25 | **4.0** | 4.0 |
| **CCGAAA** (discriminating) | 6 | 2 | 1 | 1 | 2/6 = 0.3333 | **3.0** | 2.667 |

The designed asymmetric sequence `CCGAAA` (C=2, G=1) separates the two conventions: GGF gives **3.0**, Saxonov **2.667**. The code (line 290, `expected = (C×G)/N`) yields 3.0 — confirming GGF. All other spec O/E values reproduced independently.

Worked island (400 bp CGCG): length 400 ≥ 200, GC 100% > 50%, O/E 2.0 > 0.6 → qualifies (matches M15).

### Findings / divergences
None. Description is biologically and mathematically correct, with the O/E denominator convention (C×G) matching the primary GGF source.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs`
- `FindCpGSites` lines 148–162
- `CalculateCpGObservedExpected` lines 270–292
- `FindCpGIslands` lines 297–355 (+ `CalculateGcContent` 357–358)

### Formula realised correctly?
- `FindCpGSites`: uppercases input, scans `i = 0 .. Length-2`, yields `i` where `seq[i]=='C' && seq[i+1]=='G'` — exact 5'→3' C-then-G; GpC not matched. INV-1, INV-5 hold.
- `CalculateCpGObservedExpected` (line 290–291): `expected = ((double)c * g) / Length; return expected > 0 ? cpg / expected : 0;` — exactly GGF `CpG / ((C×G)/N)`, with `expected > 0` guard covering C=0 or G=0 (INV-2, INV-3). Length<2/null guarded at line 272. **Improvement since prior report:** the `c*g` product is now computed in `double` (line 290), removing the int-overflow risk the previous report had noted as a non-defect for long windows.
- `FindCpGIslands`: sliding window, step 1, window length `minLength` (the `Math.Min` at line 314 never shortens within the loop bound `i ≤ Length-minLength`). Windows passing `gc >= minGc && cpgRatio >= minCpGRatio` are merged; the merged island is **re-validated** against GC% and O/E (lines 333–337, 348–352) and emitted only if `End-Start ≥ minLength`. Uses `>=` operators, matching Takai & Jones' "≥ 200 bp, ObsCpG/ExpCpG ≥ 0.6, %GC ≥ 50%" phrasing. INV-4 holds.

### Cross-verification table recomputed vs code
The CpG-filtered test set (102 tests incl. M9→2.0, M11→3.0, M12→4.0, M10→0.0, C1 (CCCC)→0.0, M15→(0,400,Gc=1.0,O/E=2.0)) passes. The discriminating `CCGAAA` case (GGF 3.0 vs Saxonov 2.667) was hand-traced through line 290 and matches GGF.

### Variant/delegate consistency
`FindCpGIslands` delegates O/E and GC to the same `CalculateCpGObservedExpected` / `CalculateGcContent` used standalone. `GetMethylationContext`/`FindMethylationSites` use the same C-then-G CpG rule. Consistent.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CpGDetection_Tests.cs` plus added fuzz tests (commit 55a69e3d). Assertions check exact sourced values (exact positions, `Is.EqualTo(2.0/3.0/4.0)`, exact Start/End/GcContent/CpGRatio), not tautologies; deterministic; cover all Stage-A edge cases (null, empty, single char, no C/G, GpC, adjacent, lowercase, thresholds, too-short).

### Findings / defects
None.

## Verdict & follow-ups
- Stage A: **PASS** — definition, GGF O/E formula (denominator `(C×G)/N`, confirmed distinct from Saxonov via the `CCGAAA` discriminator), thresholds, and Takai-Jones alternative all confirmed against authoritative sources.
- Stage B: **PASS** — code realises the validated formula and thresholds exactly; div-by-zero guarded; int-overflow risk now removed via `double` product; all worked examples reproduced.
- **State: CLEAN.** No defect found. No code changed.
- Build: succeeded. Tests: 102 CpG-filtered passed; full suite **18213 passed, 0 failed**.
