# Validation Report: PROBE-EVALUE-001 — Karlin-Altschul Off-Target E-value

- **Validated:** — (pending)   **Area:** MolTools
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/PROBE-EVALUE-001.md`.

## Canonical method(s)
`ComputeKarlinAltschul`, `ComputeLambdaNucleotide`

- **Source:** `src/Seqeron/Algorithms/.../ProbeDesigner.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_ProbeValidation_Tests.cs`

## Authoritative sources (to open in Stage A)
- Karlin & Altschul (1990), Altschul et al. (1990) BLAST

## Contract / invariants
R: E-value ≥ 0; M: higher bit score → lower E-value; M: larger search space → higher E-value

## Cross-check / differential oracle
- Reference: NCBI BLAST stats / published λ
- Comparison: λ≈1.374, E within tolerance

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
