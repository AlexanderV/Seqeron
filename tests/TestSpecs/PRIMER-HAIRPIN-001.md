# Test Specification: PRIMER-HAIRPIN-001

**Test Unit ID:** PRIMER-HAIRPIN-001
**Area:** MolTools
**Algorithm:** DNA Hairpin Folder + Secondary-Structure (unimolecular) Tm
**Status:** ☑ Validated — Stage A ✅ / Stage B ✅ / CLEAN (2026-06-25)
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

| # | Source | What it provides |
|---|--------|------------------|
| 1 | SantaLucia J, Hicks D (2004). Annu Rev Biophys Biomol Struct 33:415–440 | Table 1 NN stem stacks; Table 4 hairpin-loop ΔG°37 by size (ΔH°=0, ΔS°=ΔG°37·1000/310.15); Eq. 7 Jacobson-Stockmayer (coeff 2.44); Eq. 8–11 hairpin model + unimolecular Tm (no C_T term). Full PDF read this session. |
| 2 | SantaLucia J (1998). PNAS 95(4):1460–65 | Unified NN ΔH°/ΔS° (stem stacks; same values reproduced in Table 1 above). |
| 3 | primer3-py 2.3.0 `calc_hairpin` + shipped `primer3_config/{loops,triloop,tetraloop}.{dh,ds}` | Independent ntthal oracle for `CalculateHairpinThermodynamicsNtthal`; special tri/tetraloop bonus tables (triloop ±2000, tetraloop ±1100). |

## 2. Canonical Method(s)

- `PrimerDesigner.FindMostStableHairpin(string, int minStemLength=2, double loopBonusDeltaG37=0)` → `HairpinResult?`
- `PrimerDesigner.CalculateHairpinMeltingTemperature(string, int minStemLength=2, double loopBonusDeltaG37=0)` → `double`
- `PrimerDesigner.CalculateHairpinThermodynamicsNtthal(string, double sodiumMolar=0.05)` → `HairpinThermodynamics?` (bundled special tri/tetraloop bonuses)

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs` (+ `NtthalHairpin.cs`)
- **Test fixtures:** `PrimerDesigner_HairpinTm_Tests.cs`, `PrimerDesigner_HairpinSpecialLoop_Tests.cs`

## 3. Contract / Invariants

- R (range): a stable hairpin has ΔG°37 ≤ 0; if no WC stem can close a ≥3-nt loop, result is `null` / NaN Tm.
- M (monotone): a longer loop with the same stem is more destabilising (higher ΔG°37) — Jacobson-Stockmayer.
- Loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15. Stem = NN stacks only (no bimolecular init).
- Tm is unimolecular/concentration-independent: Tm = ΔH°·1000/ΔS° − 273.15 (no R·ln(C_T/x) term).
- D (deterministic): same input → same output.

## 4. Cross-check / Differential Oracle

- **Reference:** primer3-py 2.3.0 `calc_hairpin` (ntthal path) + hand-derivation from SantaLucia & Hicks 2004 Table 1/Table 4 (legacy path).
- **Comparison:** ntthal path matches primer3 to machine precision (ΔH exact; ΔS/Tm ≤1e-6). Legacy path matches hand-derivation to <1e-12.

### Worked numbers (locked in tests)
- `GGGCTTTTGCCC` (legacy Table 4): ΔH=−25.8, ΔS=−75.48486216346927, ΔG37=−2.3883700000000054, Tm=68.64038366828805 °C.
- `GGGGCGAAAGCCCC` (ntthal, GAAA tetraloop): ΔH=−40900 cal, ΔS=−114.1872884299936, ΔG37=−5484.812493437487 cal, Tm=85.03347700825856 °C (primer3 parity).
- `GGGCGAAGCCC` (ntthal, GAA triloop): ΔH=−27800 cal, Tm=84.7060915802943 °C (primer3 parity).

## 5. Validation Checklist (restored ☑)

- [x] Stage A: SantaLucia & Hicks (2004) full PDF retrieved; Table 1, Table 4, Eq. 7/11 confirmed verbatim against code constants.
- [x] Stage B: implementation reviewed; legacy path = Table 4 model, ntthal path = primer3 parity; both correct.
- [x] Independent cross-check: primer3-py 2.3.0 `calc_hairpin` (8 sequences) + hand-derivation (3 hairpins).
- [x] Coverage added: 3-bp stem, minStemLength selectivity, palindrome-no-loop, long stem+loop.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0, 0 warnings.
- [x] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
