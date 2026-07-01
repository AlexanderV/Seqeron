# Test Specification: SPLICE-MAXENT5-001

**Test Unit ID:** SPLICE-MAXENT5-001
**Area:** Splicing
**Algorithm:** MaxEntScan score5 (5' Donor)
**Status:** ☑ Complete — independently validated 2026-06-25 (Stage A ✅ / Stage B ✅ / CLEAN)
**Last Updated:** 2026-06-25

> **Validated.** Independently re-validated against Yeo & Burge (2004) and the MIT-licensed maxentpy `score5`
> reference port (see `docs/Validation/reports/SPLICE-MAXENT5-001.md`). `score5` was reimplemented in Python over
> the embedded `maxent_score5.txt` table and reproduced the published worked examples exactly
> (`cagGTAAGT`→10.858313/10.86, `gagGTAAGT`→11.078494/11.08, `taaATAAGT`→-0.116791/-0.12); the C# `ScoreDonorMaxEnt`
> matches these bit-for-bit. Full unfiltered `dotnet test Seqeron.sln` Failed: 0.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Yeo & Burge (2004), maxentpy (MIT) |

## 2. Canonical Method(s)

`ScoreDonorMaxEnt` (MaxEntScan score5)

- **Source file:** `SpliceSitePredictor.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs`

## 3. Contract / Invariants

R: score finite; D: deterministic; requires 9-nt donor window

## 4. Cross-check / Differential Oracle

- **Reference:** MaxEntScan score5.pl
- **Comparison:** exact score

## 5. Validation Checklist (restored ☑ — 2026-06-25)

- [x] Stage A: sources retrieved (Yeo & Burge 2004; maxentpy `maxent.py::score5` + `score5_matrix.txt`).
      Formula/constants (`bgd_5`, `cons1_5`, `cons2_5`, log2 factorisation, rest=fa[:3]+fa[5:]) confirmed
      against maxentpy and the published worked examples. ✅ PASS
- [x] Stage B: `ScoreDonorMaxEnt` is line-for-line equivalent to maxentpy `score5`; cross-checked vs the
      Python oracle over the embedded table — exact match. ✅ PASS
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18762 passed).
- [x] Flipped `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md` (registry row + catalog header + Status +
      Quick-Reference Completed 227→228 / Not Started 28→27). `docs/checklists/*.md` intentionally untouched.
