# Validation Report: IMMUNE-NUSVR-001 — CIBERSORT nu-SVR Immune Deconvolution (ABIS)

- **Validated:** — (pending)   **Area:** Oncology
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/IMMUNE-NUSVR-001.md`.

## Canonical method(s)
`DeconvoluteImmuneCellsNuSvr`, `LoadBundledAbisSignatureMatrix`, `LoadSignatureMatrix`

- **Source:** `src/Seqeron/Algorithms/.../ImmuneAnalyzer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ImmuneAnalyzer_ImmuneInfiltration_Tests.cs`

## Authoritative sources (to open in Stage A)
- Newman et al. (2015) CIBERSORT, Schölkopf et al. (2000) ν-SVR, Monaco et al. (2019) ABIS CC-BY

## Contract / invariants
R: fractions ≥ 0; D: deterministic; planted truth recovered

## Cross-check / differential oracle
- Reference: scikit-learn NuSVR
- Comparison: coefficients < 2e-3

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
