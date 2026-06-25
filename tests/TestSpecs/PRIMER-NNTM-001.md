# Test Specification: PRIMER-NNTM-001

**Test Unit ID:** PRIMER-NNTM-001
**Area:** MolTools
**Algorithm:** Nearest-Neighbour Salt/Mismatch/Dangling-End Tm
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
| 1 | SantaLucia (1998), Allawi & SantaLucia (1997), Owczarzy (2004/2008), Bommarito (2000) |

## 2. Canonical Method(s)

`CalculateMeltingTemperatureNN`, `CalculateMeltingTemperatureNNMismatch`

- **Source file:** `PrimerDesigner.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_NearestNeighborTm_Tests.cs`

## 3. Contract / Invariants

R: Tm finite for len ≥ 2; M: higher [Na+] → higher Tm; M: more mismatches → lower Tm

## 4. Cross-check / Differential Oracle

- **Reference:** primer3-py / Biopython MeltingTemp
- **Comparison:** Tm ±0.5°C

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
