# Test Specification: MIRNA-PCT-001

**Test Unit ID:** MIRNA-PCT-001
**Area:** MiRNA
**Algorithm:** TargetScan PCT (Branch-Length Conservation)
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

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
