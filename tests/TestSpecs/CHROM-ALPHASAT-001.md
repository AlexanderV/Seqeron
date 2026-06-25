# Test Specification: CHROM-ALPHASAT-001

**Test Unit ID:** CHROM-ALPHASAT-001
**Area:** Chromosome
**Algorithm:** Alpha-Satellite Monomer Detection
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
| 1 | Waye & Willard (1987), Henikoff et al. (2001), CENP-B box motif |

## 2. Canonical Method(s)

`DetectAlphaSatellite`, `FindCenpBBoxes`

- **Source file:** `ChromosomeAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_AlphaSatellite_Tests.cs`

## 3. Contract / Invariants

R: monomer period ≈ 171 bp; R: CENP-B boxes within monomers; D: deterministic

## 4. Cross-check / Differential Oracle

- **Reference:** known centromeric reference arrays
- **Comparison:** period + CENP-B positions agree

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
