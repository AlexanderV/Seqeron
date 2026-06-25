# Test Specification: MHC-MATRIX-001

**Test Unit ID:** MHC-MATRIX-001
**Area:** Oncology
**Algorithm:** SMM / BIMAS Matrix pMHC Prediction
**Status:** ☐ Not Started — pending independent Stage A/B re-validation
**Last Updated:** 2026-06-25

> **Stub.** This unit was added during the limitation-elimination campaign. The algorithm is implemented and
> covered by the test fixture below, but it has **not yet** been independently re-validated under the project's
> two-stage (Stage A description / Stage B implementation) protocol. This spec captures the evidence and contract
> needed to perform that validation; fill in the full TestSpec when the unit is re-validated to `☑`.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Peters et al. (2005) SMM, Parker et al. (1994) BIMAS |

## 2. Canonical Method(s)

`PredictIc50Smm`, `PredictBindingHalfLifeBimas`, `PredictAndClassifySmm`

- **Source file:** `OncologyAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs`

## 3. Contract / Invariants

R: IC50 > 0; R: BIMAS half-life ≥ 0; M: anchor match → stronger binding

## 4. Cross-check / Differential Oracle

- **Reference:** published worked examples / IEDB (caller matrix)
- **Comparison:** exact on anchor cases

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
