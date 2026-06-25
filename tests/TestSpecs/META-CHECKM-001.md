# Test Specification: META-CHECKM-001

**Test Unit ID:** META-CHECKM-001
**Area:** Metagenomics
**Algorithm:** CheckM Marker-Gene Completeness/Contamination
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
| 1 | Parks et al. (2015) CheckM, Parks et al. (2018) GTDB, Pfam CC0 |

## 2. Canonical Method(s)

`EstimateBinQualityFromMarkers`, `DetectMarkers`, `LoadBundledBacterial/ArchaealMarkerHmms`

- **Source file:** `MetagenomicsAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs`

## 3. Contract / Invariants

R: 0 ≤ completeness ≤ 100; R: contamination ≥ 0; D: deterministic

## 4. Cross-check / Differential Oracle

- **Reference:** CheckM markerSets.py
- **Comparison:** completeness/contamination exact on synthetic bin

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
