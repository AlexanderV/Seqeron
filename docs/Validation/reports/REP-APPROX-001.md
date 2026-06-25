# Validation Report: REP-APPROX-001 — Approximate (TRF) Tandem-Repeat Detection

- **Validated:** — (pending)   **Area:** Repeats
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/REP-APPROX-001.md`.

## Canonical method(s)
`FindApproximateTandemRepeats`, `ComputeBernoulliStatistics`

- **Source:** `src/Seqeron/Algorithms/.../RepeatFinder.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_ApproximateTandemRepeats_Tests.cs`

## Authoritative sources (to open in Stage A)
- Benson (1999) TRF

## Contract / invariants
R: percent-matches ∈ [0,100]; R: score ≥ Minscore (50); D: deterministic

## Cross-check / differential oracle
- Reference: TRF (Benson) on benchmark repeats
- Comparison: consensus + match/indel% agree

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
