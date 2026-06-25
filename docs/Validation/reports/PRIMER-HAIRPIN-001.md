# Validation Report: PRIMER-HAIRPIN-001 — DNA Hairpin Folder + Secondary-Structure Tm

- **Validated:** 2026-06-25   **Area:** MolTools
- **Canonical method(s):** `PrimerDesigner.FindMostStableHairpin`, `PrimerDesigner.CalculateHairpinMeltingTemperature`, `PrimerDesigner.CalculateHairpinThermodynamicsNtthal` (bundled special tri/tetraloop bonuses)
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
- `FindMostStableHairpin(string, int minStemLength=2, double loopBonusDeltaG37=0)` → `HairpinResult?`
- `CalculateHairpinMeltingTemperature(string, int minStemLength=2, double loopBonusDeltaG37=0)` → `double`
- `CalculateHairpinThermodynamicsNtthal(string, double sodiumMolar=0.05)` → `HairpinThermodynamics?` (full Primer3 ntthal monomer DP; auto-applies bundled triloop/tetraloop bonus tables)

Source: `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs` (lines ~1183–1368, ~1659–1712) and `NtthalHairpin.cs`.
Tests: `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_HairpinTm_Tests.cs`, `PrimerDesigner_HairpinSpecialLoop_Tests.cs`.

## Authoritative sources opened THIS session
1. **SantaLucia J, Hicks D (2004). _The thermodynamics of DNA structural motifs._ Annu Rev Biophys Biomol Struct 33:415–440** — full PDF retrieved and text-extracted (Duke mirror). Confirmed:
   - **Table 1** (p.420) Watson-Crick NN ΔH°/ΔS° (1 M NaCl) — verbatim match to `NnUnifiedParams`; Initiation +0.2/−5.7; Terminal A·T +2.2/+6.9; Symmetry 0.0/−1.4.
   - **Table 4** (p.427) hairpin-loop ΔG°37 increments by size — verbatim match to `HairpinLoopInitiationDeltaG` (3→3.5, 4→3.5, 5→3.3, 6→4.0, 7→4.2, 8→4.3, 9→4.5, 10→4.6, 12→5.0, 14→5.1, 16→5.3, 18→5.5, 20→5.7, 25→6.1, 30→6.3).
   - Table 4 footnotes: "All loop ΔH° parameters are assumed to equal zero"; "ΔS° = ΔG°37 × 1000/310.15"; loops <3 sterically prohibited; loops of length 3/4 require special triloop/tetraloop corrections (supplementary).
   - **Eq. 7** Jacobson-Stockmayer extrapolation `ΔG°37(loop-n) = ΔG°37(loop-x) + 2.44·R·310.15·ln(n/x)`, coefficient **2.44** (preferred over older 1.75) — verbatim match.
   - **Eq. 8–10** (p.428) hairpin total: stem NN stacks (Table 1) + loop energy; length-3 adds a +0.5 closing-A·T penalty + triloop bonus; length-4 adds triloop/tetraloop + terminal mismatch; length ≥4 adds terminal mismatch.
   - **Eq. 11** (p.428) unimolecular two-state hairpin Tm = ΔH°·1000/ΔS° − 273.15 (**no** strand-concentration term) — verbatim match.
2. **primer3-py 2.3.0** (`pip install primer3-py==2.3.0`) — `calc_hairpin` used as the independent ntthal oracle; its shipped `primer3_config/{loops,triloop,tetraloop}.{dh,ds}` tables inspected directly (triloop CGAAG/GGAAC = −2000 dh / 0 ds; tetraloop CGAAAG/GGGGAC = −1100 dh / 0 ds; loops.ds hairpin column size3 = −10.31, size4 = −11.6).

## Stage A — Description

- **Source quality** — PASS. Primary literature (the actual 2004 Annu Rev paper) opened and read; every constant the code claims to source from it is present verbatim. primer3-py config tables inspected at the file level.
- **Formula correctness** — PASS. Two distinct, both-correct models are implemented and clearly delineated:
  - Legacy `FindMostStableHairpin` / `CalculateHairpinMeltingTemperature` = SantaLucia & Hicks (2004) Table 4 ΔG°37 hairpin-loop model + Table 1 stem stacks, **unimolecular** Tm (Eq. 11), bimolecular duplex-initiation term correctly **excluded**.
  - `CalculateHairpinThermodynamicsNtthal` = full Primer3 ntthal monomer DP (reproduces `calc_hairpin`) with the bundled sequence-specific triloop/tetraloop bonus tables auto-applied.
- **Definitions/conventions** — PASS. 5'→3'; loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15 (destabilising → ΔS°<0, sign correct); min loop = 3; ΔG°37 = ΔH° − 310.15·ΔS°/1000; Tm has no C_T term. The 5'-arm-dinucleotide NN stacking is exact because the NN table is reverse-complement symmetric (verified for all 16 steps).
- **Edge-case semantics** — PASS. No complementary stem → no hairpin (poly-A); loop <3 prohibited; non-ACGT/empty/null → none; palindrome with 0-nt loop → none. All sourced.
- **Independent cross-check** — see numbers below.
- **Invariants** — PASS. ΔG°37 ≤ 0 for the canonical stable cases; longer loop → higher (less negative) ΔG°37 (Jacobson-Stockmayer); deterministic.

### Documented model boundary (not a defect)
The legacy folder bundles the stem-stack + Table 4 loop-initiation **core** only. The length-3 closing-A·T penalty (+0.5), the length-3/4 special triloop/tetraloop bonus, and the ≥4 terminal-mismatch increment (Eq. 8–10 supplementary terms) are **not** baked into the legacy folder; they are exposed as an opt-in additive `loopBonusDeltaG37`. This is explicitly disclosed in the source ("NOT bundled (honest residual)") and is consistent with the source — the core values it does compute are exact. The **bundled** versions of those special-loop terms are provided by the separate `CalculateHairpinThermodynamicsNtthal` path, validated below against primer3.

## Stage B — Implementation

- **Code path reviewed** — `PrimerDesigner.cs` `HairpinLoopInitiationDeltaG` / `HairpinLoopDeltaG` (Jacobson-Stockmayer) / `FindMostStableHairpin` / `CalculateHairpinMeltingTemperature` / `CalculateHairpinThermodynamicsNtthal`; `NtthalHairpin.cs`.
- **Formula realised correctly** — PASS. Stem stacks summed over 5'-arm dinucleotides (exact by NN symmetry), no bimolecular init term, loop ΔH°=0, loop ΔS°=−ΔG·1000/310.15, ΔG°37 and unimolecular Tm exactly as Eq. 10/11. Greedy inward stem extension from each WC outermost pair, MFE over all candidates — a single perfect WC stem (no internal mismatch/bulge), which matches the model scope.
- **Edge cases in code** — PASS. `minStemLength<2`, empty/null, non-ACGT → null; loop<3 skipped; no-stem → null; `loopBonusDeltaG37` added to loop ΔG with matching ΔS contribution.
- **Variant/delegate consistency** — PASS. `CalculateHairpinMeltingTemperature` delegates to `FindMostStableHairpin` and applies Eq. 11. ntthal path returns native cal/mol converted to the library kcal/mol convention.
- **Numerical robustness** — PASS. No overflow/div-by-zero on stated ranges; loop-31..∞ via Jacobson-Stockmayer.
- **Tests are real** — PASS, strengthened. Expected values trace to source/hand-derivation/primer3, not code echoes; tight tolerances (1e-9 / 1e-7); special-loop values asserted to primer3-py 2.3.0 to machine precision. No skips, no weakened assertions.

### Cross-verification — primer3-py 2.3.0 `calc_hairpin` (mv=50, dv=0, dntp=0, dna_conc=50 nM) vs C# `CalculateHairpinThermodynamicsNtthal`

| Seq | Loop | primer3 ΔH (cal) | primer3 ΔS | primer3 ΔG37 (cal) | primer3 Tm °C | C# match |
|-----|------|----------------:|-----------:|-------------------:|--------------:|:--------:|
| GGGGCGAAAGCCCC | GAAA tetraloop | −40900 | −114.1872884299936 | −5484.812493437487 | 85.03347700825856 | ✅ |
| GGGGGGGACCCCC | GGGA tetraloop | −34000 | −94.1872884299836 | −4787.81249344059 | 87.8328944728006 | ✅ |
| GGGCGAAGCCC | GAA triloop | −27800 | −77.68485895331574 | −3706.040995629125 | 84.7060915802943 | ✅ |
| GGGGGAACCCC | GAA triloop | −26000 | −73.18485895331571 | −3301.715995629134 | 82.11474153055735 | ✅ |
| GGGCTTTTGCCC | TTTT (non-special) | −32400 | −94.58485895332572 | −3064.5059956260297 | 69.39954078842845 | ✅ |
| GGGCAAAAAGCCC | 5-nt (non-special) | −30100 | −87.74485895332572 | −2885.9319956260306 | 69.89004085311882 | ✅ |
| AAAAAAAAAAAA | — | structure_found=False | | | | ✅ (null) |
| GCGC | — | structure_found=False | | | | ✅ (null) |

Every value matches to machine precision (asserted exactly for ΔH, ParityTol=1e-6 for ΔS/Tm).

### Cross-verification — hand-derivation (legacy Table 4 model)

Canonical `GGGCTTTTGCCC` (4-bp GGGC/GCCC stem, 4-nt TTTT loop), SantaLucia & Hicks 2004 Table 1 + Table 4:
- Stem ΔH° = GG+GG+GC = −8.0−8.0−9.8 = **−25.8** kcal/mol; Stem ΔS° = −19.9−19.9−24.4 = **−64.2** e.u.
- Loop-4 ΔG°37 = 3.5 → loop ΔS° = −3.5·1000/310.15 = **−11.28486216346929** e.u.
- Total ΔH° = **−25.8**; Total ΔS° = **−75.48486216346927**; ΔG°37 = −25.8 − 310.15·(−75.48486…)/1000 = **−2.3883700000000054**.
- **Tm = −25.8·1000/−75.48486216346927 − 273.15 = 68.64038366828805 °C** (Eq. 11, no C_T term).

Hand values reproduced independently this session and matched by the code to <1e-12. New tests also lock a 3-bp-stem hairpin (`GCCAAAGGC`: ΔH=−17.8, ΔS=−55.58486216346929, ΔG37=−0.5603550000000013, Tm=47.08107204353689 °C) and a 6-bp-stem/10-nt-loop hairpin (ΔH=−40.0, ΔS=−114.33153312913106, ΔG37=−4.540075000000002), all hand-derived from Table 1/Table 4.

### Tool divergence (documented)
primer3 `ntthal` uses the SantaLucia (1998) NN set + its own loop model (loops.ds hairpin column, e.g. size4 = −11.6 e.u.) **plus** terminal-mismatch (`tstack2`) and special tri/tetraloop bonuses — so `calc_hairpin` on `GGGCTTTTGCCC` returns ΔH=−32400 cal, Tm=69.40 °C, whereas the **legacy SantaLucia & Hicks 2004 Table 4** model returns ΔH=−25.8 kcal, Tm=68.64 °C. These are two different (both correctly sourced) models, not an error: the C# `CalculateHairpinThermodynamicsNtthal` matches primer3 exactly, while `FindMostStableHairpin` matches the 2004 Annu-Rev Table 4 hand-derivation exactly. The library exposes both and documents the boundary.

## Test-quality audit
- `PrimerDesigner_HairpinTm_Tests.cs` — 14 tests (10 original + 4 added this session: 3-bp stem, minStemLength selectivity, palindrome-no-loop, long stem+loop). Covers thermodynamics, unimolecular Tm, concentration-independence, no-hairpin, loop-size table selection, min-loop violation, invalid input, minStem<2 guard, Jacobson-Stockmayer, opt-in bonus, competing stems, varying stem/loop length, palindrome, very long.
- `PrimerDesigner_HairpinSpecialLoop_Tests.cs` — special tri/tetraloop bonus path (2 tetraloop + 2 triloop), non-special-loop regression (4-nt, 5-nt), homopolymer/too-short/invalid → null, legacy-model-unchanged regression. All primer3-sourced.
- No green-washing found: tolerances are tight, every expected number is sourced/hand-derived/primer3, all Stage-A edge cases covered.

## Verdict & follow-ups
- **Stage A: ✅ PASS** — description and constants match SantaLucia & Hicks (2004) verbatim; both models correctly sourced; model boundary honestly documented.
- **Stage B: ✅ PASS** — code realises both formulas exactly; cross-checks match primer3-py 2.3.0 to machine precision and hand-derivation to <1e-12.
- **State: ✅ CLEAN.** No defect. Test suite strengthened (+4 hand-derived coverage tests for varying stem length, minStemLength selectivity, palindrome, and long structures). Full unfiltered `dotnet test Seqeron.sln -c Debug`: Failed 0, 0 warnings.
- Intermolecular self-/hetero-dimer Tm is a **separate** unit (PRIMER-DIMER-001) and out of scope here.
