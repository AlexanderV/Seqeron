# Validation Report: CHROM-HOR-001 — Higher-Order Repeat (HOR) Detection

- **Validated:** — (pending)   **Area:** Chromosome
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/CHROM-HOR-001.md`.

## Canonical method(s)
`DetectHigherOrderRepeat`

- **Source:** `src/Seqeron/Algorithms/.../ChromosomeAnalyzer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs`

## Authoritative sources (to open in Stage A)
- McNulty & Sullivan (2018), Alkan et al. (2007)

## Contract / invariants
R: inter-HOR identity ≥ intra-monomer identity; R: HOR period = k×monomer; D: deterministic

## Cross-check / differential oracle
- Reference: known HOR arrays (D-region)
- Comparison: period + copy number agree

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
