# Validation Report: MIRNA-CLASSIFY-001 — Pre-miRNA Structure-Feature Classifier

- **Validated:** — (pending)   **Area:** MiRNA
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/MIRNA-CLASSIFY-001.md`.

## Canonical method(s)
`ClassifyPreMiRna` (logistic over MFE/AMFE/MFEI/GC/%paired)

- **Source:** `src/Seqeron/Algorithms/.../MiRnaAnalyzer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs`

## Authoritative sources (to open in Stage A)
- Bonnet et al. (2004), miRBase (public domain), Zhang (2006) MFEI

## Contract / invariants
R: probability ∈ [0,1]; D: deterministic; threshold split positive/negative

## Cross-check / differential oracle
- Reference: held-out miRBase vs shuffled (AUC)
- Comparison: AUC ≈ 1.0 on held-out set

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
