# Test Specification: PROTMOTIF-HMM-001

**Test Unit ID:** PROTMOTIF-HMM-001
**Area:** ProteinMotif
**Algorithm:** Plan7 Profile-HMM Domain Search
**Status:** ☐ Not Started — pending independent Stage A/B re-validation
**Last Updated:** 2026-06-25

> **Stub.** This unit was added during the limitation-elimination campaign. The algorithm is implemented and
> covered by the test fixture below, but it has **not yet** been independently re-validated under the project's
> two-stage (Stage A description / Stage B implementation) protocol. This spec captures the evidence and contract
> needed to perform that validation; fill in the full TestSpec when the unit is re-validated to `☑`.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Eddy (1998, 2011) HMMER/Plan7, Pfam |

## 2. Canonical Method(s)

`Plan7ProfileHmm.Viterbi/Forward`, `FindDomainsByHmm`, `FindDomainEnvelopes`

- **Source file:** `Plan7ProfileHmm.cs / ProteinMotifFinder.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`

## 3. Contract / Invariants

R: Forward ≥ Viterbi (log-odds); R: E-value ≥ 0; D: deterministic given profile

## 4. Cross-check / Differential Oracle

- **Reference:** pyhmmer / hmmsearch
- **Comparison:** bit score ±1e-3, same envelopes

## 5. Validation Checklist (to restore ☑)

- [ ] Stage A: retrieve every source above; confirm formula/constants against the publication's worked example.
- [ ] Stage B: review the implementation against the source; cross-check vs the reference oracle.
- [ ] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [ ] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.
