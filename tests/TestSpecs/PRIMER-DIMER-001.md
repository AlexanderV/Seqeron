# Test Specification: PRIMER-DIMER-001

**Test Unit ID:** PRIMER-DIMER-001
**Area:** MolTools
**Algorithm:** ntthal Self/Hetero-Dimer Tm (thermodynamic alignment)
**Status:** ☑ Complete — independently validated 2026-06-25; re-reviewed 2026-10-01 (campaign 2026-09, B07: Stage A ⚠ → fixed, Stage B ❌ → fixed, see `docs/Validation/review-2026-09/B07.md` F11)
**Last Updated:** 2026-10-01

---

## 1. Evidence Summary

| # | Source | What it provides |
|---|--------|------------------|
| 1 | SantaLucia J, Hicks D (2004) *Annu Rev Biophys* 33:415–440 | Unified NN ΔH°/ΔS° (Table 1); bimolecular Tm Eq. 3 (x=1 self-comp / x=4 non); 0.368 salt correction (Eq. 5) |
| 2 | Untergasser A et al. (2012) *Nucleic Acids Res* 40:e115 | Primer3/ntthal thermodynamic alignment for bimolecular duplexes |
| 3 | primer3 `thal.c` + `primer3_config/*.dh,*.ds` | DP (`fillMatrix`/`LSH`/`RSH`/`maxTM`/`calc_bulge_internal`/`traceback`/`calcDimer`) + parameter tables |
| 4 | primer3-py 2.3.0 `calc_homodimer`/`calc_heterodimer` | Reference ΔH/ΔS/ΔG/Tm oracle (mv=50, dv=0, dntp=0, dna=50 nM) |
| 5 | primer3-py 2.3.1 `thal.c` (tag v2.3.1) + `calc_heterodimer`/`calc_homodimer`/`calc_end_stability` (all arguments, `output_structure=True`) | Bit-faithful oracle for ANY/END1/END2, mv/dv/dntp/dna_conc/temp_c/max_loop, ASCII structure, length limits |

## 2. Canonical Method(s)

`FindMostStableDimer`, `CalculateDimerMeltingTemperature`, `CalculateSelfDimerMeltingTemperature`,
`CalculateDimerThermodynamicsNtthal` (full DP; overloads with alignment mode / Mg²⁺ / dNTP / temp_c / max_loop),
`CalculateDimerStructureNtthal` (+ ASCII structure), internal `NtthalDimer.Run`.

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`, `NtthalDimer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_DimerTm_Tests.cs`

## 3. Contract / Invariants

- **R (range):** `null` / NaN only when ntthal finds no structure; ntthal's optimum may be a single pair
  with negative Tm and/or positive ΔG (e.g. primer3-py `calc_heterodimer('A','T')`: Tm −437.07 °C,
  ΔG +2087.8 cal/mol; END2 `CGGGCG`/`CGTCCGAATCAGTTGTGAT`: ΔG +1967.85 cal/mol) — reported, as primer3-py does.
- **L (limits):** both strands > 60 nt or either > 10 000 nt → `ArgumentException` (thal.c CHECK_ERROR);
  max_loop outside 0–30 → `ArgumentOutOfRangeException` (primer3-py setter). Opt-in overloads with a trailing
  `maxAlignLength` (60–10 000, else `ArgumentOutOfRangeException`) replace 60 by that THAL_MAX_ALIGN — identical to
  thal.c compiled with `-DTHAL_MAX_ALIGN=…` (audit round 3, A3-9, F55; `PrimerDesigner_NtthalMaxAlign_Tests`).
- **T (temperature):** temp_c changes only ΔG (= ΔH − (temp_c+273.15)·ΔS), not Tm / ΔH / ΔS.
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

### 4.1 Re-validation 2026-10-01 (PRIMER-DIMER-001, primer3-py 2.3.1)

- Root cause of the ≈1.4 % residual mismatch (ΔG off by ±13.8…1014 cal/mol, e.g. 208.5 cal/mol):
  thal.c `LSH`/`RSH` keep `T1 = −∞` unless a dangling-end branch is entered, so the bare A·T
  candidate wins over the tstack2 candidate; the port recomputed T1 from the tstack2 candidate.
- After the fix: 8000 random pairs (5–60 nt; random, self-complementary, GC-rich, homopolymer runs,
  36–60 nt; ANY/END1/END2; default and random mv/dv/dntp/dna_conc/temp_c/max_loop): 0 mismatches,
  max |ΔTm| = 0, max |ΔG| = 7.3e-12 cal/mol (before: 30/4000 mismatches at the defaults, max |ΔG| 1014).
  ASCII structure: 1034/1034 ANY pairs identical to `ascii_structure_lines`.
- Tests: `CalculateDimerThermodynamicsNtthal_FormerDiscrepancies_MatchPrimer3Py` (6),
  `_TemperatureAndMaxLoop_MatchPrimer3Py`, `CalculateDimerStructureNtthal_AsciiStructure_MatchesPrimer3Py`,
  `CalculateDimerThermodynamicsNtthal_ThalLimits`.

## 5. Validation Checklist (restored ☑)

- [x] Stage A: every source retrieved; formula/constants confirmed against SantaLucia & Hicks Table 1 + thal.c.
- [x] Stage B: implementation reviewed against thal.c; cross-checked vs primer3-py 2.3.0 oracle.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Genomics.Tests 18741 passed).
- [x] Flipped `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
