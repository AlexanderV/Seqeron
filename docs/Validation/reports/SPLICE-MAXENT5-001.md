# Validation Report: SPLICE-MAXENT5-001 — MaxEntScan score5 (5' Donor)

- **Validated:** — (pending)   **Area:** Splicing
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/SPLICE-MAXENT5-001.md`.

## Canonical method(s)
`ScoreDonorMaxEnt` (MaxEntScan score5)

- **Source:** `src/Seqeron/Algorithms/.../SpliceSitePredictor.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs`

## Authoritative sources (to open in Stage A)
- Yeo & Burge (2004), maxentpy (MIT)

## Contract / invariants
R: score finite; D: deterministic; requires 9-nt donor window

## Cross-check / differential oracle
- Reference: MaxEntScan score5.pl
- Comparison: exact score

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
