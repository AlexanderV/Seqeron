# Validation Report: PRIMER-HAIRPIN-001 — DNA Hairpin Folder + Secondary-Structure Tm

- **Validated:** — (pending)   **Area:** MolTools
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/PRIMER-HAIRPIN-001.md`.

## Canonical method(s)
`FindMostStableHairpin`, `CalculateHairpinMeltingTemperature`

- **Source:** `src/Seqeron/Algorithms/.../PrimerDesigner.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_HairpinTm_Tests.cs`

## Authoritative sources (to open in Stage A)
- SantaLucia (1998), SantaLucia & Hicks (2004), UNAFold DNA params

## Contract / invariants
R: ΔG of best hairpin ≤ 0 (or none found); M: longer stem → more negative ΔG; D: deterministic

## Cross-check / differential oracle
- Reference: mfold / UNAFold (if obtainable)
- Comparison: ΔG ±0.2 kcal/mol, Tm ±1°C

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
