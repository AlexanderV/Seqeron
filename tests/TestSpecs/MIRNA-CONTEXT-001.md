# Test Specification: MIRNA-CONTEXT-001

**Test Unit ID:** MIRNA-CONTEXT-001
**Area:** MiRNA
**Algorithm:** TargetScan context++ Scoring
**Status:** ☑ Validated — Stage A ✅ PASS / Stage B ✅ PASS / ✅ CLEAN (2026-06-25)
**Last Updated:** 2026-06-25

> Independently re-validated 2026-06-25 under the two-stage protocol; see
> `docs/Validation/reports/MIRNA-CONTEXT-001.md`. All Agarwal_2015_parameters.txt coefficients verified
> byte-exact; every feature computation/scaling confirmed against `targetscan_70_context_scores.pl`
> (`getAgarwalContribution`, `getLocalAU_contribution`, `get_sRNA1_8_contributions`,
> `getSite8_contribution`, `getSA_contribution`, `get3primePairingContribution`, `getMinDist`,
> `get_len3UTR`, `getOffset6mer`); PCT sigmoid against `targetscan_70_BL_PCT.pl` `calculatePCTthisBL`.
> SA derives from the Turner-2004 McCaskill partition function (Z_open/Z), not MFE. Full-transcript
> features (TA/SPS/Len_ORF/ORF8m) and PCT are a documented caller-supplied boundary (computed verbatim
> when supplied), not a limitation.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Agarwal et al. (2015) TargetScan context++ |

## 2. Canonical Method(s)

`ScoreTargetSiteContextPlusPlus` (+ SA accessibility wiring)

- **Source file:** `MiRnaAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs`

## 3. Contract / Invariants

R: context++ score ≤ 0 (more negative = stronger); D: deterministic

## 4. Cross-check / Differential Oracle

- **Reference:** targetscan_70_context_scores.pl
- **Comparison:** computable subset byte-exact

## 5. Validation Checklist (to restore ☑)

- [x] Stage A: retrieved Agarwal 2015 eLife + `Agarwal_2015_parameters.txt` + `targetscan_70_context_scores.pl` + `targetscan_70_BL_PCT.pl`; confirmed coefficients/formula/scaling byte-exact and hand-computed the locally-computable subset (intercept + Local_AU + sRNA1/8 + Site8 + SA + 3P + Min_dist + Len_3UTR + Off6m + PCT).
- [x] Stage B: reviewed the implementation against the perl reference; cross-checked the computable subset byte-exact (CTX-001/004/011 + all PCT cases recomputed independently).
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Seqeron.Genomics.Tests 18762 passed).
- [x] Flipped `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md`.
