# Test Specification: SPLICE-MAXENT5-001

**Test Unit ID:** SPLICE-MAXENT5-001
**Area:** Splicing
**Algorithm:** MaxEntScan score5 (5' Donor)
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

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
