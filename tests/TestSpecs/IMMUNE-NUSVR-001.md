# Test Specification: IMMUNE-NUSVR-001

**Test Unit ID:** IMMUNE-NUSVR-001
**Area:** Oncology
**Algorithm:** CIBERSORT nu-SVR Immune Deconvolution (ABIS)
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
| 1 | Newman et al. (2015) CIBERSORT, Schölkopf et al. (2000) ν-SVR, Monaco et al. (2019) ABIS CC-BY |

## 2. Canonical Method(s)

`DeconvoluteImmuneCellsNuSvr`, `LoadBundledAbisSignatureMatrix`, `LoadSignatureMatrix`

- **Source file:** `ImmuneAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ImmuneAnalyzer_ImmuneInfiltration_Tests.cs`

## 3. Contract / Invariants

R: fractions ≥ 0; D: deterministic; planted truth recovered

## 4. Cross-check / Differential Oracle

- **Reference:** scikit-learn NuSVR
- **Comparison:** coefficients < 2e-3

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
