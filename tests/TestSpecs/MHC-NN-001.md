# Test Specification: MHC-NN-001

**Test Unit ID:** MHC-NN-001
**Area:** Oncology
**Algorithm:** MHCflurry Pan-Allele NN Binding Affinity
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
| 1 | O'Donnell et al. (2018, 2020) MHCflurry, Apache-2.0 |

## 2. Canonical Method(s)

`MhcflurryAffinityPredictor.PredictIc50`, ensemble geometric-mean combiner

- **Source file:** `MhcflurryAffinityPredictor.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MhcflurryAffinityPredictor_PredictIc50_Tests.cs`

## 3. Contract / Invariants

R: IC50 > 0; R: ensemble within member range; D: deterministic given weights

## 4. Cross-check / Differential Oracle

- **Reference:** mhcflurry 2.1.5 (models_class1_pan)
- **Comparison:** IC50 < 0.03%

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
