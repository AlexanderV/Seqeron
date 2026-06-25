# Validation Report: PRIMER-NNTM-001 — Nearest-Neighbour Salt/Mismatch/Dangling-End Tm

- **Validated:** — (pending)   **Area:** MolTools
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/PRIMER-NNTM-001.md`.

## Canonical method(s)
`CalculateMeltingTemperatureNN`, `CalculateMeltingTemperatureNNMismatch`

- **Source:** `src/Seqeron/Algorithms/.../PrimerDesigner.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_NearestNeighborTm_Tests.cs`

## Authoritative sources (to open in Stage A)
- SantaLucia (1998), Allawi & SantaLucia (1997), Owczarzy (2004/2008), Bommarito (2000)

## Contract / invariants
R: Tm finite for len ≥ 2; M: higher [Na+] → higher Tm; M: more mismatches → lower Tm

## Cross-check / differential oracle
- Reference: primer3-py / Biopython MeltingTemp
- Comparison: Tm ±0.5°C

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
