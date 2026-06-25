# Test Specification: PROBE-EVALUE-001

**Test Unit ID:** PROBE-EVALUE-001
**Area:** MolTools
**Algorithm:** Karlin-Altschul Off-Target E-value
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
| 1 | Karlin & Altschul (1990), Altschul et al. (1990) BLAST |

## 2. Canonical Method(s)

`ComputeKarlinAltschul`, `ComputeLambdaNucleotide`

- **Source file:** `ProbeDesigner.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs`

## 3. Contract / Invariants

R: E-value ≥ 0; M: higher bit score → lower E-value; M: larger search space → higher E-value

## 4. Cross-check / Differential Oracle

- **Reference:** NCBI BLAST stats / published λ
- **Comparison:** λ≈1.374, E within tolerance

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
