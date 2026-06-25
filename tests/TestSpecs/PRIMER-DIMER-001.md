# Test Specification: PRIMER-DIMER-001

**Test Unit ID:** PRIMER-DIMER-001
**Area:** MolTools
**Algorithm:** ntthal Self/Hetero-Dimer Tm (thermodynamic alignment)
**Status:** ☑ Complete — independently validated 2026-06-25 (Stage A PASS / Stage B PASS / CLEAN)
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

| # | Source | What it provides |
|---|--------|------------------|
| 1 | SantaLucia J, Hicks D (2004) *Annu Rev Biophys* 33:415–440 | Unified NN ΔH°/ΔS° (Table 1); bimolecular Tm Eq. 3 (x=1 self-comp / x=4 non); 0.368 salt correction (Eq. 5) |
| 2 | Untergasser A et al. (2012) *Nucleic Acids Res* 40:e115 | Primer3/ntthal thermodynamic alignment for bimolecular duplexes |
| 3 | primer3 `thal.c` + `primer3_config/*.dh,*.ds` | DP (`fillMatrix`/`LSH`/`RSH`/`maxTM`/`calc_bulge_internal`/`traceback`/`calcDimer`) + parameter tables |
| 4 | primer3-py 2.3.0 `calc_homodimer`/`calc_heterodimer` | Reference ΔH/ΔS/ΔG/Tm oracle (mv=50, dv=0, dntp=0, dna=50 nM) |

## 2. Canonical Method(s)

`FindMostStableDimer`, `CalculateDimerMeltingTemperature`, `CalculateSelfDimerMeltingTemperature`,
`CalculateDimerThermodynamicsNtthal` (full DP), internal `NtthalDimer.Run`.

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`, `NtthalDimer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_DimerTm_Tests.cs`

## 3. Contract / Invariants

- **R (range):** a dimer's optimal structure has ΔG ≤ 0, else no structure (`null` / NaN).
- **M (monotonicity):** a longer/stronger complementary core has the lower (more stable) ΔG / higher Tm.
- **D (determinism):** identical inputs ⇒ identical outputs.
- **INV:** self-dimer(S) ≡ hetero-dimer(S, S). x=1 only when both strands are RC palindromes.
- **Salt/conc:** lower [Na⁺] lowers Tm (0.368·N·ln[Na⁺]); higher C_T raises bimolecular Tm.

## 4. Cross-check / Differential Oracle

- **Reference:** primer3-py 2.3.0 `calc_homodimer`/`calc_heterodimer`.
- **Result:** C# reproduces primer3-py to machine precision on contiguous-WC optima **and** on
  internal-loop (2×2, 3×3, mixed), single-base-bulge, and terminal-overhang optima; hand-derived
  SantaLucia & Hicks values match to 1e-9. See `docs/Validation/reports/PRIMER-DIMER-001.md` for the
  full per-case table.

## 5. Validation Checklist (restored ☑)

- [x] Stage A: every source retrieved; formula/constants confirmed against SantaLucia & Hicks Table 1 + thal.c.
- [x] Stage B: implementation reviewed against thal.c; cross-checked vs primer3-py 2.3.0 oracle.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Genomics.Tests 18741 passed).
- [x] Flipped `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
