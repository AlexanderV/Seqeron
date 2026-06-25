# Validation Report: RNA-ACCESS-001 — McCaskill Unpaired (Accessibility) Probabilities

- **Validated:** — (pending)   **Area:** RnaStructure
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/RNA-ACCESS-001.md`.

## Canonical method(s)
`CalculateUnpairedProbabilities`, `CalculateRegionUnpairedProbability`

- **Source:** `src/Seqeron/Algorithms/.../RnaSecondaryStructure.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_UnpairedProbabilities_Tests.cs`

## Authoritative sources (to open in Stage A)
- McCaskill (1990), RNAplfold/Bernhart (2006)

## Contract / invariants
R: 0 ≤ P_unpaired ≤ 1; M: longer region → lower P_unpaired; D: deterministic

## Cross-check / differential oracle
- Reference: brute-force ensemble enumeration (small n)
- Comparison: equal P_unpaired ±1e-9

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
