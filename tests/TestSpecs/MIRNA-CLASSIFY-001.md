# Test Specification: MIRNA-CLASSIFY-001

**Test Unit ID:** MIRNA-CLASSIFY-001
**Area:** MiRNA
**Algorithm:** Pre-miRNA Structure-Feature Classifier
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

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
