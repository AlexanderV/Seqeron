# Validation Report: META-CHECKM-001 — CheckM Marker-Gene Completeness/Contamination

- **Validated:** — (pending)   **Area:** Metagenomics
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/META-CHECKM-001.md`.

## Canonical method(s)
`EstimateBinQualityFromMarkers`, `DetectMarkers`, `LoadBundledBacterial/ArchaealMarkerHmms`

- **Source:** `src/Seqeron/Algorithms/.../MetagenomicsAnalyzer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs`

## Authoritative sources (to open in Stage A)
- Parks et al. (2015) CheckM, Parks et al. (2018) GTDB, Pfam CC0

## Contract / invariants
R: 0 ≤ completeness ≤ 100; R: contamination ≥ 0; D: deterministic

## Cross-check / differential oracle
- Reference: CheckM markerSets.py
- Comparison: completeness/contamination exact on synthetic bin

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
