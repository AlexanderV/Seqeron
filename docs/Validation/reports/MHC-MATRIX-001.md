# Validation Report: MHC-MATRIX-001 — SMM / BIMAS Matrix pMHC Prediction

- **Validated:** — (pending)   **Area:** Oncology
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/MHC-MATRIX-001.md`.

## Canonical method(s)
`PredictIc50Smm`, `PredictBindingHalfLifeBimas`, `PredictAndClassifySmm`

- **Source:** `src/Seqeron/Algorithms/.../OncologyAnalyzer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_ClassifyMhcBinding_Tests.cs`

## Authoritative sources (to open in Stage A)
- Peters et al. (2005) SMM, Parker et al. (1994) BIMAS

## Contract / invariants
R: IC50 > 0; R: BIMAS half-life ≥ 0; M: anchor match → stronger binding

## Cross-check / differential oracle
- Reference: published worked examples / IEDB (caller matrix)
- Comparison: exact on anchor cases

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
