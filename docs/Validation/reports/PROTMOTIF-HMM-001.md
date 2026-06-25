# Validation Report: PROTMOTIF-HMM-001 — Plan7 Profile-HMM Domain Search

- **Validated:** — (pending)   **Area:** ProteinMotif
- **Stage A verdict:** ⬜ pending
- **Stage B verdict:** ⬜ pending
- **State:** ⬜ pending first validation

> **Stub.** Net-new algorithm added during the limitation-elimination campaign. The code is implemented and
> covered by the fixture below, but this unit has **not yet** been independently validated under
> [VALIDATION_PROTOCOL.md](../VALIDATION_PROTOCOL.md). This report is a placeholder to be completed when the
> unit is validated; see also `tests/TestSpecs/PROTMOTIF-HMM-001.md`.

## Canonical method(s)
`Plan7ProfileHmm.Viterbi/Forward`, `FindDomainsByHmm`, `FindDomainEnvelopes`

- **Source:** `src/Seqeron/Algorithms/.../Plan7ProfileHmm.cs / ProteinMotifFinder.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`

## Authoritative sources (to open in Stage A)
- Eddy (1998, 2011) HMMER/Plan7, Pfam

## Contract / invariants
R: Forward ≥ Viterbi (log-odds); R: E-value ≥ 0; D: deterministic given profile

## Cross-check / differential oracle
- Reference: pyhmmer / hmmsearch
- Comparison: bit score ±1e-3, same envelopes

## Stage A — Description
_Pending._

## Stage B — Implementation
_Pending._

## Verdict
⬜ **Pending first validation.**
