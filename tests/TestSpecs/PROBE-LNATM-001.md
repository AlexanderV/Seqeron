# Test Specification: PROBE-LNATM-001

**Test Unit ID:** PROBE-LNATM-001
**Area:** MolTools
**Algorithm:** LNA-Adjusted NN Tm + MGB Probe Design
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
| 1 | McTigue (2004) LNA NN, Kutyavin (2000) MGB |

## 2. Canonical Method(s)

`CalculateMeltingTemperatureNNLna`, `EvaluateMgbProbeDesign`

- **Source file:** `ProbeDesigner.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_LnaTm_Tests.cs`

## 3. Contract / Invariants

R: each LNA substitution does not lower Tm; D: deterministic; MGB rules return boolean+reasons

## 4. Cross-check / Differential Oracle

- **Reference:** MELTING 5
- **Comparison:** Tm ±0.2°C

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
