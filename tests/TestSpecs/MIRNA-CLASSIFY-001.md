# Test Specification: MIRNA-CLASSIFY-001

**Test Unit ID:** MIRNA-CLASSIFY-001
**Area:** MiRNA
**Algorithm:** Pre-miRNA Structure-Feature Classifier
**Status:** ☑ Complete — independently validated 2026-06-25 (Stage A ✅ / Stage B ✅ / CLEAN)
**Last Updated:** 2026-06-25

> **Validated.** Independently re-validated under the two-stage protocol; see
> `docs/Validation/reports/MIRNA-CLASSIFY-001.md`. miRBase positives verified against the live database, the
> AMFE/MFEI/GC/%paired features and the logistic score reproduced by hand to full double precision, di-shuffle
> negatives (two independent seeds) score below threshold, and the held-out AUC = 1.0. Hardening tests CL13/CL14
> added.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Bonnet et al. (2004), miRBase (public domain), Zhang (2006) MFEI |

## 2. Canonical Method(s)

`ClassifyPreMiRna` (logistic over MFE/AMFE/MFEI/GC/%paired)

- **Source file:** `MiRnaAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs`

## 3. Contract / Invariants

R: probability ∈ [0,1]; D: deterministic; threshold split positive/negative

## 4. Cross-check / Differential Oracle

- **Reference:** held-out miRBase vs shuffled (AUC)
- **Comparison:** AUC ≈ 1.0 on held-out set

## 5. Validation Checklist (to restore ☑)

- [x] Stage A: retrieved Bonnet 2004 (di-shuffle null model), Zhang 2006 (AMFE/MFEI), miRBase (sequences verified live); formulas confirmed.
- [x] Stage B: implementation matches the source; features and logistic score reproduced by hand to full precision; held-out AUC = 1.0.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Genomics 18771 passed, 0 warnings).
- [x] Flipped `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md` (registry + catalog + Quick-Reference). (docs/checklists/* intentionally untouched.)
