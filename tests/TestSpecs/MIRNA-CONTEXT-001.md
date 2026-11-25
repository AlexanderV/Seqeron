# Test Specification: MIRNA-CONTEXT-001

**Test Unit ID:** MIRNA-CONTEXT-001
**Area:** MiRNA
**Algorithm:** TargetScan context++ Scoring
**Status:** ☑ Validated — Stage A 🟡 PASS-WITH-NOTES / Stage B FAIL→FIXED (B13 re-review 2026-11-25, F15–F18)
**Last Updated:** 2026-11-25

> Independently re-validated 2026-06-25 under the two-stage protocol; see
> `docs/Validation/reports/MIRNA-CONTEXT-001.md`. All Agarwal_2015_parameters.txt coefficients verified
> byte-exact; every feature computation/scaling confirmed against `targetscan_70_context_scores.pl`
> (`getAgarwalContribution`, `getLocalAU_contribution`, `get_sRNA1_8_contributions`,
> `getSite8_contribution`, `getSA_contribution`, `get3primePairingContribution`, `getMinDist`,
> `get_len3UTR`, `getOffset6mer`); PCT sigmoid against `targetscan_70_BL_PCT.pl` `calculatePCTthisBL`.
> SA derives from the Turner-2004 McCaskill partition function (Z_open/Z), not MFE. Full-transcript
> features (TA/SPS/Len_ORF/ORF8m) and PCT are a documented caller-supplied boundary (computed verbatim
> when supplied), not a limitation. **Update 2026-11-25 (B13 F33/F34):** TA/SPS default to the bundled
> TargetScan `TA_SPS_by_seed_region.txt`; Len_ORF/ORF8m are derived from `ContextPlusPlusInputs.OrfSequence`
> (`targetscan_count_8mers.pl`). Caller inputs left: PCT alignment evidence (tree + species with site) and
> the exact RNAplfold SA value.

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

- R: `ContextScore` (TargetScan's reported value) ≤ the site-type ceiling (8mer −0.03, 7mer-m8 −0.02, 7mer-A1 −0.01, 6mer 0); `null` when the 1-based site start is < 15 (`too_close`).
- R: `ContextScore` = round3(Σ round3(contribution)), capped; `ContextScorePartial` = unrounded, uncapped sum.
- R: Local_AU fraction is rounded to 3 decimals before min-max scaling; Off6m counts non-overlapping occurrences.
- R: SA uses the caller-supplied RNAplfold value when given (then equal to TargetScan); 0 contribution and omitted when the probability is 0 or the 14-nt window has no `_lunp` row; a site within 7 nt of the 3' end reads the last row (`grep -A7 | tail -1`).
- R: TA / SPS = the bundled TargetScan row of the seed region (nt 2–8; SPS column 2 for 8mer/7mer-m8, column 3 for 7mer-A1/6mer) unless the caller supplies them; omitted only when the seed region is not 7 nt over ACGU.
- R: Len_ORF / ORF8m from `OrfSequence`: gaps `-`/`.` removed, length = remaining characters; ORF8m = case-insensitive, non-overlapping count of `revcomp(seed)+"A"` (T as U); caller `OrfLength` / `Orf8mCount` win.
- D: deterministic.

## 4. Cross-check / Differential Oracle

- **Reference:** targetscan_70_context_scores.pl (TargetScan 7.0, nsoranzo/targetscan main), run with perl 5.38.2; RNAplfold values from ViennaRNA 2.8.0 `pfl_fold_up(seq, 20, 80, 40)`.
- **Comparison:** `MiRnaAnalyzer_ContextPlusPlusTargetScanParity_Tests` — 181-nt UTR, let-7a, sites 8mer 31–38 / 6mer 82–87 / 7mer-A1 125–131 from `targetscan_70.pl`: all 14 printed contribution columns and the context++ score (−0.742 / −0.118 / −0.257) equal the perl output when the RNAplfold value is supplied.
- **Kit oracle (F33/F34):** `TargetScan7_context_scores.zip` — `targetscan_70_context_scores_output.txt` (675 rows): all TA/SPS contributions from the bundled table; ORF lengths + 8mer counts of 5 kit ORF rows (14 family seeds) equal `ORF_Sequences_sample.lengths.txt` / `ORF_8mer_counts_sample.txt`; Len_ORF/ORF8m of the 91 output rows of those genes. Worked example with only `OrfSequence` + RNAplfold value: all 14 columns + scores equal the perl run fed by `targetscan_count_8mers.pl`; a second ORF (2535 nt, 3 8mers) equals the perl's ORF columns and scores (−0.929 / −0.231 / −0.358).
- **Known reference defects not reproduced:** Off6m is always 1 in the public perl (match against a hash reference); Site8 for 6mer reads the base two positions 5' of the one opposite miRNA nt 8 (22-nt subsequence padded with two N's). Without a supplied RNAplfold value SA is an approximation (one 80-nt window, no 40-nt span limit).

## 5. Validation Checklist (to restore ☑)

- [x] Stage A: retrieved Agarwal 2015 eLife + `Agarwal_2015_parameters.txt` + `targetscan_70_context_scores.pl` + `targetscan_70_BL_PCT.pl`; confirmed coefficients/formula/scaling byte-exact and hand-computed the locally-computable subset (intercept + Local_AU + sRNA1/8 + Site8 + SA + 3P + Min_dist + Len_3UTR + Off6m + PCT).
- [x] Stage B: reviewed the implementation against the perl reference; cross-checked the computable subset byte-exact (CTX-001/004/011 + all PCT cases recomputed independently).
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Seqeron.Genomics.Tests 18762 passed).
- [x] Flipped `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md`.
- [x] B13 re-review 2026-11-25: reference perl run (above); F15 reported score (rounding, ceiling, too_close), F16 Local_AU rounding, F17 SA (supplied RNAplfold value, zero / 3'-end rules; estimate documented as approximation), F18 Off6m non-overlapping. Tests: `MiRnaAnalyzer_ContextPlusPlusTargetScanParity_Tests` (22).
- [x] B13 WP-4/5 2026-11-25: F33 TA/SPS bundled (`TryGetTargetScanTaSps`, default when not supplied), F34 Len_ORF/ORF8m from `OrfSequence` (`CountTargetScanOrf8merSites`). Tests: `MiRnaAnalyzer_ContextPlusPlusTargetScanParity_Tests` 22 → 45 cases (+ kit data `MiRnaContextTargetScanKitData.cs`); hand-derived CTX-001..004 / CTX-007 and the Algebraic/Fuzzing worked-example pins re-locked with the bundled TA/SPS terms.
