# Test Specification: PRIMER-HAIRPIN-001

**Test Unit ID:** PRIMER-HAIRPIN-001
**Area:** MolTools
**Algorithm:** DNA Hairpin Folder + Secondary-Structure (unimolecular) Tm
**Status:** ☑ Validated — Stage A ✅ / Stage B ✅ / CLEAN (2026-06-25); re-reviewed 2026-10-01 (review-2026-09 B07): ntthal engine made bit-exact to primer3-py 2.3.1 (F13–F14)
**Last Updated:** 2026-10-01

---

## 1. Evidence Summary

| # | Source | What it provides |
|---|--------|------------------|
| 1 | SantaLucia J, Hicks D (2004). Annu Rev Biophys Biomol Struct 33:415–440 | Table 1 NN stem stacks; Table 4 hairpin-loop ΔG°37 by size (ΔH°=0, ΔS°=ΔG°37·1000/310.15); Eq. 7 Jacobson-Stockmayer (coeff 2.44); Eq. 8–11 hairpin model + unimolecular Tm (no C_T term). Full PDF read this session. |
| 2 | SantaLucia J (1998). PNAS 95(4):1460–65 | Unified NN ΔH°/ΔS° (stem stacks; same values reproduced in Table 1 above). |
| 3 | primer3-py 2.3.0 `calc_hairpin` + shipped `primer3_config/{loops,triloop,tetraloop}.{dh,ds}` | Independent ntthal oracle for `CalculateHairpinThermodynamicsNtthal`; special tri/tetraloop bonus tables (triloop ±2000, tetraloop ±1100). |
| 4 | primer3-py **2.3.1** vendored `primer3/src/libprimer3/thal.c` + `thal.h` + `primer3/thermoanalysis.pyx` (raw.githubusercontent.com/libnano/primer3-py/v2.3.1/…) | The type-4 path ported line by line: `thal`, `initMatrix2`, `fillMatrix2`, `maxTM2`, `CBI`, `calc_bulge_internal2`, `calc_hairpin`, `RSH`, `Ss`/`Hs`, `calc_terminal_bp`, `END5_1..4`, `max5`, `tracebacku`, `equal`, `calcHairpin`, `drawHairpin`, `THAL_MAX_ALIGN`, `temp_c`/`max_loop`. |

## 2. Canonical Method(s)

- `PrimerDesigner.FindMostStableHairpin(string, int minStemLength=2, double loopBonusDeltaG37=0)` → `HairpinResult?`
- `PrimerDesigner.CalculateHairpinMeltingTemperature(string, int minStemLength=2, double loopBonusDeltaG37=0)` → `double`
- `PrimerDesigner.CalculateHairpinThermodynamicsNtthal(string, double sodiumMolar=0.05)` → `HairpinThermodynamics?` (bundled special tri/tetraloop bonuses; dv = dntp = 0)
- `PrimerDesigner.CalculateHairpinThermodynamicsNtthal(string, mv, dv, dntp[, temperatureCelsius, maxLoop])` → `HairpinThermodynamics?` (= primer3-py `calc_hairpin`)
- `PrimerDesigner.CalculateHairpinStructureNtthal(string, mv=0.05, dv=0.0015, dntp=0.0006, temperatureCelsius=37, maxLoop=30)` → `NtthalHairpinStructure?` (+ `ascii_structure_lines`)
- `PrimerDesigner.FindMostStableHairpin(string, StructureModel, double? mv, dv, dntp, temperatureCelsius, int? maxLoop)` / `CalculateHairpinMeltingTemperature(string, StructureModel, …)` (audit round 3, A3-13 / F57): `SingleHelix` = the default single-stem core (no condition arguments), `Ntthal` = `calc_hairpin` mapped onto `HairpinResult` (StemStart = 5′-most paired base, StemEnd = its partner, StemLength = pairs of that stem across bulges/internal loops, LoopSize = unpaired bases inside its innermost pair)

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs` (+ `NtthalHairpin.cs`)
- **Test fixtures:** `PrimerDesigner_HairpinTm_Tests.cs`, `PrimerDesigner_HairpinSpecialLoop_Tests.cs`

## 3. Contract / Invariants

- R (range): a stable hairpin has ΔG°37 ≤ 0; if no WC stem can close a ≥3-nt loop, result is `null` / NaN Tm.
- M (monotone): a longer loop with the same stem is more destabilising (higher ΔG°37) — Jacobson-Stockmayer.
- Loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15. Stem = NN stacks only (no bimolecular init).
- Tm is unimolecular/concentration-independent: Tm = ΔH°·1000/ΔS° − 273.15 (no R·ln(C_T/x) term).
- D (deterministic): same input → same output.
- ntthal path: identical to primer3-py 2.3.1 `calc_hairpin` (Tm, ΔG at `temp_c`, ΔH, ΔS, ASCII structure, `structure_found`); > 60 nt → `ArgumentException` (opt-in `maxAlignLength` overloads: > maxAlignLength, 60–10 000; = thal.c compiled with `-DTHAL_MAX_ALIGN=…`, 61–120-nt values locked in `PrimerDesigner_NtthalMaxAlign_Tests`, A3-9 / F55); `maxLoop` outside 0–30 → `ArgumentOutOfRangeException`; non-ACGT / empty → `null`.

## 4. Cross-check / Differential Oracle

- **Reference:** primer3-py 2.3.0 `calc_hairpin` (ntthal path) + hand-derivation from SantaLucia & Hicks 2004 Table 1/Table 4 (legacy path).
- **Comparison:** ntthal path matches primer3 to machine precision (ΔH exact; ΔS/Tm ≤1e-6). Legacy path matches hand-derivation to <1e-12.

- **2026-10-01 (B07):** 9000 random oligos (5–60 nt; random, palindromic, GC-rich, homopolymer runs, designed hairpins with mismatches/bulges) at the defaults and at random mv/dv/dntp/temp_c/max_loop — 0 mismatches (≤ 1e−6) in Tm/ΔG/ΔH/ΔS and ASCII structure (pre-fix port: 185/1000 and 575/3000 mismatches). `DesignPrimers` vs `design_primers`: 1800/1800 random templates identical (was 1733/1800).

### Worked numbers (locked in tests)
- `GGGCTTTTGCCC` (legacy Table 4): ΔH=−25.8, ΔS=−75.48486216346927, ΔG37=−2.3883700000000054, Tm=68.64038366828805 °C.
- `GGGGCGAAAGCCCC` (ntthal, GAAA tetraloop): ΔH=−40900 cal, ΔS=−114.1872884299936, ΔG37=−5484.812493437487 cal, Tm=85.03347700825856 °C (primer3 parity).
- `GGGCGAAGCCC` (ntthal, GAA triloop): ΔH=−27800 cal, Tm=84.7060915802943 °C (primer3 parity).
- Former discrepancies (calc_hairpin defaults): `GGGAGACAGTAGTCGCCCAT` Tm 64.43690682436392 (old 69.31), `TGTTGAATATCAGCG` ΔG +530.6133132321556 cal, `TTTGCCACTAATAATATGATCAACCGGAGGGTCTCCATT` Tm 83.58095561848427 (old 139.59), a 53-mer with no structure (old: 32.64 °C) — `PrimerDesigner_HairpinTm_Tests.CalculateHairpinThermodynamicsNtthal_FormerDiscrepancies_MatchPrimer3Py`, `_FormerFalseStructure_IsNoStructure`, `_ConditionsTemperatureAndMaxLoop_MatchPrimer3Py`, `CalculateHairpinStructureNtthal_AsciiStructure_MatchesPrimer3Py`, `_ThalLimits`.

### Audit round 3, A3-13 (F57)
- Single-stem core: every stem length closing a ≥ 3-nt loop is scored (INV-1/INV-2; was maximal extension only):
  `GGGGCCCC` → 0/6, 2 bp, loop 3, ΔH −8.0, ΔS −31.18486216346929, ΔG°37 +1.6719850000000012, Tm −16.61527625921667 °C
  (was null); `GCGCGCGCGC` → 0/9, 3 bp, loop 4, ΔH −20.4, ΔS −62.884862163469286, ΔG°37 −0.8962600000000016,
  Tm 51.25239666853645 °C (was null; primer3-py `calc_hairpin` draws the same 3-bp/4-nt shape). Hand-derived from
  Table 1 + Table 4. Table 4 loop values re-checked against primer3 `loops.ds` (hairpin column, −ΔS·310.15/1000:
  3.4985/3.4985/3.300/3.998/4.199/4.299/4.497/4.600/5.000/5.099/5.297/5.499/5.697/6.098/6.299 for 3…30) and the
  Jacobson–Stockmayer fill (n = 11: 4.7433 vs 4.742; n = 13: 5.120 vs 5.117).
- `StructureModel.Ntthal` vs `calc_hairpin(…, output_structure=True)`: 8000 random oligos (seeds 1–2; 5–60 nt; 40 % random,
  60 % designed stem-loops with mismatches/indels and flanks; half at the defaults, half at random
  mv/dv/dntp/temp_c/max_loop) — 0 differences (ΔH/ΔS/ΔG/Tm |Δ| = 0; span/stem/loop = `ascii_structure_lines`).
  Locked: `GGGCTTTTGCCC` 0/11/4/4, −32.4, −93.71131562841765, −3.335435457846266, Tm 72.59266493570397;
  `CCCTGAGTCCGAGGAGAGGGT` 0/19/6/4 (internal loop), ΔG +0.40853085538899177, Tm 34.24603584502296;
  `TTTGCCACTAATAATATGATCAACCGGAGGGTCTCCATT` 21/31/3/3, Tm 83.58095561848427 — `PrimerDesigner_StructureModel_Tests`.

## 5. Validation Checklist (restored ☑)

- [x] Stage A: SantaLucia & Hicks (2004) full PDF retrieved; Table 1, Table 4, Eq. 7/11 confirmed verbatim against code constants.
- [x] Stage B: implementation reviewed; legacy path = Table 4 model, ntthal path = primer3 parity; both correct.
- [x] Independent cross-check: primer3-py 2.3.0 `calc_hairpin` (8 sequences) + hand-derivation (3 hairpins).
- [x] Coverage added: 3-bp stem, minStemLength selectivity, palindrome-no-loop, long stem+loop.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0, 0 warnings.
- [x] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
