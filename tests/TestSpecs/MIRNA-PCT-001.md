# Test Specification: MIRNA-PCT-001

**Test Unit ID:** MIRNA-PCT-001
**Area:** MiRNA
**Algorithm:** TargetScan PCT (Branch-Length Conservation)
**Status:** ☑ Complete — independently validated (Stage A ✅ / Stage B ✅ / CLEAN), 2026-06-25
**Last Updated:** 2026-06-25

> **Validated.** Independently re-validated under the two-stage protocol against Friedman et al. (2009) Genome Res
> 19:92 and the TargetScan reference Perl `targetscan_70_BL_PCT.pl` (`calculatePCTthisBL` logistic + `getBranchLength`
> + `PCT_parameters` tables, all retrieved verbatim). See `docs/Validation/reports/MIRNA-PCT-001.md`. No code defect;
> the PCT fixture was strengthened (49→56 cases) to cover all four site-type parameter branches, the monotonicity
> invariant, the untruncated Bls=0 logistic floor, and a real published-parameter cross-check.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Friedman et al. (2009) PCT, TargetScan |

## 2. Canonical Method(s)

PCT (branch-length-score → logistic) from caller alignment+tree

- **Source file:** `MiRnaAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs`

## 3. Contract / Invariants

R: PCT ∈ [0,1]; M: higher branch length → higher PCT; D: deterministic

## 4. Cross-check / Differential Oracle

- **Reference:** Friedman 2009 logistic worked example
- **Comparison:** PCT within tolerance

## 5. Validation Checklist (restored ☑)

- [x] Stage A: retrieved `targetscan_70_BL_PCT.pl` + `8mer_PCT_parameters.txt` + Friedman 2009; formula
      `PCT=b0+b1/(1+e^(−b2·Bls+b3))` (truncated at 0) and BLS (minimal connecting-subtree branch length) confirmed verbatim.
- [x] Stage B: implementation matches the source; cross-checked Bls {A,B}=3.0/{A,C}=7.0/{A,B,C,D}=12.0/{A}=0.0 and
      PCT(3)=0.952574…/PCT(7)=0.999088…/PCT(0)=0.5/truncation/real miR-30 8mer row vs an independent reference.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18769 passed); PCT fixture 56 passed.
- [x] Flipped `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md` (registry row + catalog header + Quick-Reference counts).
      (`docs/checklists/*.md` intentionally NOT modified per the validation-campaign deliverable rules.)
