# Test Specification: REP-APPROX-001

**Test Unit ID:** REP-APPROX-001
**Area:** Repeats
**Algorithm:** Approximate (TRF) Tandem-Repeat Detection
**Status:** ☑ Complete — independently validated 2026-06-25 (Stage A ✅ / Stage B ✅ / CLEAN)
**Last Updated:** 2026-06-25

> **Validated.** Independently re-validated under the two-stage protocol on 2026-06-25 against Benson (1999) and
> the official TRF docs (definitions + parameters), retrieved that session. Scoring (+2/−7/−7), majority-rule
> consensus, %matches/%indels "between adjacent copies", the Bernoulli (PM/PI) model, and Minscore=50 all match
> the source. Cross-checked by hand-computation on controlled perfect / substitution / indel tracts (TRF binary
> not installable). No defect; coverage strengthened with 3 edge-case tests. See
> `docs/Validation/reports/REP-APPROX-001.md`.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Benson (1999) TRF |

## 2. Canonical Method(s)

`FindApproximateTandemRepeats`, `ComputeBernoulliStatistics`

- **Source file:** `RepeatFinder.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_ApproximateTandemRepeats_Tests.cs`

## 3. Contract / Invariants

R: percent-matches ∈ [0,100]; R: score ≥ Minscore (50); D: deterministic

## 4. Cross-check / Differential Oracle

- **Reference:** TRF (Benson) on benchmark repeats
- **Comparison:** consensus + match/indel% agree

## 5. Validation Checklist (restored ☑)

- [x] Stage A: sources retrieved (Benson 1999 PMC148217; TRF definitions/parameters pages); weights +2/−7/−7, majority-rule consensus, %matches/%indels "between adjacent copies", PM=.80/PI=.10, Minscore confirmed.
- [x] Stage B: implementation reviewed vs source; hand-computed cross-check (perfect/substitution/indel + Bernoulli) reproduced exactly by the code. TRF binary not installable → hand-computation is the oracle (protocol-permitted).
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Seqeron.Genomics.Tests 18778 passed).
- [x] Flipped `☐ → ☑` in root `ALGORITHMS_CHECKLIST_V2.md` (registry + catalog header + Quick-Reference counts). (`docs/checklists/*.md` intentionally NOT modified per session scope.)
