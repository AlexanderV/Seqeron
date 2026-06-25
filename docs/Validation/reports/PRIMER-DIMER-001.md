# Validation Report: PRIMER-DIMER-001 — ntthal Self/Hetero-Dimer Tm

- **Validated:** — (pending)   **Area:** MolTools
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/PRIMER-DIMER-001.md`.

## Canonical method(s)
`FindMostStableDimer`, `CalculateDimerMeltingTemperature`, `CalculateSelfDimerMeltingTemperature`

- **Source:** `src/Seqeron/Algorithms/.../PrimerDesigner.cs / NtthalDimer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_DimerTm_Tests.cs`

## Authoritative sources (to open in Stage A)
- SantaLucia & Hicks (2004), primer3 ntthal

## Contract / invariants
R: dimer ΔG ≤ 0 or none; M: longer complementary run → lower ΔG; D: deterministic

## Cross-check / differential oracle
- Reference: primer3-py calc_homodimer/heterodimer
- Comparison: ΔG/Tm to machine precision on contiguous optima

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
