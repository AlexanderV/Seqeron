# Test Specification: PRIMER-NNTM-001

**Test Unit ID:** PRIMER-NNTM-001
**Area:** MolTools
**Algorithm:** Nearest-Neighbour Salt/Mismatch/Dangling-End Tm
**Status:** ☑ Validated (Stage A ✅ / Stage B ✅ / CLEAN — 2026-06-25)
**Last Updated:** 2026-06-25

> Independently re-validated under the two-stage protocol on 2026-06-25; see
> `docs/Validation/reports/PRIMER-NNTM-001.md`. CLEAN, no defect.

---

## 1. Evidence Summary

| # | Source | Used for |
|---|--------|----------|
| 1 | SantaLucia & Hicks (2004) *Annu Rev Biophys* 33:415 — Table 1, Eq. 3, Eq. 5 | Unified NN ΔH°/ΔS° (the constants the code uses), Tm equation, entropy salt correction |
| 2 | SantaLucia (1998) *PNAS* 95:1460 | Original unified set (init convention differs from 2004; code uses 2004) |
| 3 | Allawi & SantaLucia (1997/1998); Peyret et al. (1999) | Internal single-mismatch NN (= Biopython DNA_IMM1) |
| 4 | Bommarito, Peyret & SantaLucia (2000) *NAR* 28:1929 | Single dangling-end NN (= Biopython DNA_DE1) |
| 5 | Owczarzy et al. (2004) *Biochemistry* 43:3537 | Monovalent Na⁺ 1/Tm correction (Biopython method 6) |
| 6 | Owczarzy et al. (2008) *Biochemistry* 47:5336 | Divalent Mg²⁺/dNTP correction (Biopython method 7) |
| 7 | primer3-py 2.3.0 `calc_tm`; Biopython 1.85 `Tm_NN` | Differential oracles |

## 2. Canonical Method(s)

`CalculateMeltingTemperatureNN`, `CalculateMeltingTemperatureNNMismatch`
(+ `CalculateNearestNeighborThermodynamics`, `CalculateNearestNeighborThermodynamicsMismatch`)

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_NearestNeighborTm_Tests.cs` (25 tests)

## 3. Contract / Invariants

- R: Tm finite for ≥ 2 ACGT bases; NaN/null for empty/null/non-ACGT/len<2/unequal-length/tandem-mismatch.
- Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15, R = 1.9872, x = 4 (non-self) / 1 (self-complementary).
- M: higher [Na⁺] → higher Tm (Owczarzy 2004 monotonic); adding Mg²⁺ → higher Tm (Owczarzy 2008).
- M: an internal mismatch lowers Tm vs the repaired perfect duplex.
- ID: a fully-paired duplex through the mismatch path = the perfect-match path exactly.

## 4. Cross-check / Differential Oracle (exact numbers, 2026-06-25)

| Case | C# | Biopython `Tm_NN` (DNA_NN4) | primer3 | Reconciliation |
|---|---|---|---|---|
| ATGCATGC no-salt | 30.4338 | 30.4389 | — | R 1.9872 vs 1.987 |
| ATGCATGC Owczarzy2004 50 mM | 18.1900 | 18.1947 (m6) | 18.6227 | primer3 uses 1998 table |
| EcoRI CGCGAATTCGCG no-salt | 61.1452 | — | 61.2532 (Na=1M) | within ±0.5 °C |
| CGCGAATTCGCG Owczarzy2008 Na50/Mg3 mM | 55.4498 | 55.4529 (m7) | — | R constant only |
| MM1 internal G·T Tm | −6.4061 | −6.3997 | — | R constant only |

NN/IMM/DE tables verified verbatim vs Biopython DNA_NN4 / DNA_IMM1 / DNA_DE1 (16/50/32 entries — all match).

## 5. Validation Checklist (restored ☑)

- [x] Stage A: all sources retrieved; Tm equation, NN/init/terminal/symmetry, IMM, DE, and both Owczarzy
      corrections confirmed against the publications and reproduced by hand.
- [x] Stage B: implementation reviewed; tables verbatim vs Biopython; values match oracles within the
      reconciled R-constant / table-version differences.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed 0 (Seqeron.Genomics.Tests 18737 passed).
- [x] Flipped `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the `docs/checklists/*.md`.
