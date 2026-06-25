# Validation Report: PROBE-LNATM-001 — LNA-Adjusted NN Tm + MGB Probe Design

- **Validated:** — (pending)   **Area:** MolTools
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/PROBE-LNATM-001.md`.

## Canonical method(s)
`CalculateMeltingTemperatureNNLna`, `EvaluateMgbProbeDesign`

- **Source:** `src/Seqeron/Algorithms/.../ProbeDesigner.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_LnaTm_Tests.cs`

## Authoritative sources (to open in Stage A)
- McTigue (2004) LNA NN, Kutyavin (2000) MGB

## Contract / invariants
R: each LNA substitution does not lower Tm; D: deterministic; MGB rules return boolean+reasons

## Cross-check / differential oracle
- Reference: MELTING 5
- Comparison: Tm ±0.2°C

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
