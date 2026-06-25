# Validation Report: MHC-NN-001 — MHCflurry Pan-Allele NN Binding Affinity

- **Validated:** — (pending)   **Area:** Oncology
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/MHC-NN-001.md`.

## Canonical method(s)
`MhcflurryAffinityPredictor.PredictIc50`, ensemble geometric-mean combiner

- **Source:** `src/Seqeron/Algorithms/.../MhcflurryAffinityPredictor.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/MhcflurryAffinityPredictor_PredictIc50_Tests.cs`

## Authoritative sources (to open in Stage A)
- O'Donnell et al. (2018, 2020) MHCflurry, Apache-2.0

## Contract / invariants
R: IC50 > 0; R: ensemble within member range; D: deterministic given weights

## Cross-check / differential oracle
- Reference: mhcflurry 2.1.5 (models_class1_pan)
- Comparison: IC50 < 0.03%

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
