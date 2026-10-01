# Test Specification: REP-APPROX-001

**Test Unit ID:** REP-APPROX-001
**Area:** Repeats
**Algorithm:** Approximate (TRF) Tandem-Repeat Detection + TRF Bernoulli statistics
**Status:** ☑ Complete — re-validated 2026-09-30 (campaign 2026-09, batch B04; Stage A corrected, Stage B fixed); TRF parameters / outputs added 2026-10-01 (B04 audit WP6); TRF detection pipeline completed 2026-10-01 (B04 audit WP7)
**Last Updated:** 2026-10-01 (WP7)

> **2026-09 review.** The 2026-06 validation had no TRF binary and hand-derived its expectations from a
> window-vs-consensus model that is not TRF's. TRF 4.10.0 was compiled from source and used as the oracle:
> % matches / % indels are between ADJACENT copies, copy number counts aligned consensus columns, the alignment
> is a local wraparound DP, N never matches, and overlapping repeats of different periods are reported. See
> `docs/Evidence/REP-APPROX-001-Evidence.md`.

---

## 1. Evidence Summary

| # | Source | Used for |
|---|--------|----------|
| 1 | Benson G (1999) NAR 27:573 (content via TRF README) | model, WDP, sum-of-heads R(d,k,PM), tuple sizes, consensus, period definition, Minscore 50, PM .80 / PI .10 |
| 2 | TRF 4.10.0 README (github.com/Benson-Genomics-Lab/TRF) | parameters, table / alignment explanation, redundancy, MaxPeriod ≤ 2000, test_seqs expected tables |
| 3 | TRF 4.10.0 source, compiled (`trf seq.fa 2 7 7 80 10 <min> <maxp> -h -d`) + instrumented copy | numeric oracle for every expected value |
| 4 | TRF 4.10.0 README parameters `-m` / `-f` / `-r` / `-l` / `-ngs`, PM/PI data note; compiled TRF with non-default `Match Mismatch Delta PM PI Minscore MaxPeriod`, `-m`, `-f`, `-r`; a TRF build taking `-l` in bp (WP6) | parameter sets, masked file, flanks, alignment rows, redundancy-off, `-l` cap |
| 5 | TRF 4.10.0 README "Apparent Size Distribution" (definition of S, conditioning on the sum-of-heads criterion, 95 %, example PM .75 / k 5 / d 100 → 56), "Narrow Band Alignment" (band radius Δd_max, recentring), "Multiple Reporting …", What's New 4.04 (wider forward band) / 4.07b (alignment continues through zero); TRF source read for behaviour only (`new_meet_criteria_3`, `search_for_range_in_bestperiodlist`, `narrowbandwrap`, `get_narrowband_pair_alignment_with_copynumber`, `waitdata80/75` used as oracle) (WP7) | detection criteria, best-period list, band WDP |

## 2. Canonical Method(s)

| Method | Type |
|---|---|
| `FindApproximateTandemRepeats(DnaSequence, minPeriod = 1, maxPeriod = 6, minScore = 50)` | Canonical |
| `FindApproximateTandemRepeats(string, …)` | Overload (case-insensitive; null/empty → empty) |
| `ComputeBernoulliStatistics(string repeatTract, int period, double expectedMatchProbability = 0.80)` | Canonical |
| `FindApproximateTandemRepeats(string | DnaSequence, TandemRepeatsFinderParameters, minPeriod = 1)` | Canonical (full TRF parameter set, WP6) |
| `MaskApproximateTandemRepeats(string, TandemRepeatsFinderParameters? = null, softMask = false)` / `(string, IEnumerable<ApproximateTandemRepeatResult>, softMask)` | TRF `-m` masked sequence (+ soft mask) |
| `ApproximateTandemRepeatResult.EntropyTrf / AlignedSequence / AlignedConsensus / LeftFlank / RightFlank` | TRF entropy column, alignment rows, `-f` flanks |
| `TrfSumOfHeadsCriterion(int d)`, `TrfSumOfHeadsCriterion(int d, int pm)` (internal) | Helper (tested; PM 80 and 75) |
| `TrfApparentSize(int d, int pm)`, `TrfApparentSizeOffset(int d, int pm)` (internal) | Helper (apparent-size criterion y and its window offset max(d,20) − y − 1; tested, WP7) |

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`
- **Test fixtures:** `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_ApproximateTandemRepeats_Tests.cs`,
  `Unit/Analysis/RepeatFinder_TrfParameters_Tests.cs` (WP6: P1–P14), `Unit/Analysis/RepeatFinder_TrfDetection_Tests.cs` (WP7: D1–D9)
- **Other tiers:** `Fuzzing/RepeatApproxFuzzTests.cs`, `Properties/RepeatFinderProperties.cs` (REP-APPROX-001 region),
  `Metamorphic/RepeatsMetamorphicTests.cs`, `Combinatorial/RepeatsCombinatorialTests.cs`

## 3. Contract / Invariants

| ID | Invariant |
|---|---|
| INV-1 | `AlignmentScore ≥ minScore`; `CopyNumber ≥ 1.9` (consensus ≤ 50; ramp to 1.8; 1.8 above 100) |
| INV-2 | `PercentMatches + PercentIndels ≤ 100`, each ∈ [0,100]; composition sum ≤ 100; `Entropy` ∈ [0,2] |
| INV-3 | `minPeriod ≤ Period ≤ maxPeriod`; `|Consensus| = ConsensusSize ≥ 1` (may differ from `Period`) |
| INV-4 | No redundant pair remains (≥ 90 % overlap and same period with no higher score, or multiple period with score ≤ 1.1×); ordered by `Start` |
| INV-5 | Deterministic; case-insensitive; `DnaSequence` and string overloads agree |
| INV-6 | N / non-ACGT never match; all-N input → no repeat |
| INV-7 | Bernoulli: `Matches + Mismatches + Indels = BernoulliTrials`; PM, PI ∈ [0,1]; equal to TRF's adjacent-copy counts |
| INV-8 | `EntropyTrf` = −Σ p_b log₂ p_b over A/C/G/T with p_b = count_b / region length (N counted in the length); `EntropyTrf = Entropy` when the region is pure A/C/G/T |
| INV-9 | `AlignedSequence` and `AlignedConsensus` have equal length; `AlignedSequence` without '-' = the region; `AlignedConsensus` without '-' starts with `Consensus` |
| INV-10 | Mask: output length = input length; exactly the positions of reported repeats become `N` (soft: lower case), all others unchanged |
| INV-11 | Flanks: `LeftFlank` = the min(FlankLength, Start) symbols before the repeat, `RightFlank` = up to FlankLength symbols after it; null when FlankLength = 0 |
| INV-12 | Apparent-size offset: 1 ≤ `TrfApparentSizeOffset(d, pm)` ≤ max(d, 20) − 1 for every d = 1..2000, PM 80 / 75 |
| V-2 | `TandemRepeatsFinderParameters.Validate`: weights ≥ 1, PM ∈ {75, 80}, PI 1..100, MinScore ≥ 1, MaxPeriod 1..2000, MaxRepeatLength ≥ 1, FlankLength ≥ 0 → else `ArgumentOutOfRangeException`; null parameters / sequence → `ArgumentNullException`; mask with a repeat outside the sequence → `ArgumentOutOfRangeException` |
| V-1 | Eager `ArgumentOutOfRangeException`: `minPeriod < 1`, `maxPeriod < minPeriod`, `maxPeriod > 2000`, `minScore < 1` (both overloads, also for empty input); `ArgumentNullException` for null `DnaSequence` / tract; Bernoulli: `period ∉ 1..2000`, PM ∉ [0,1] or NaN, tract < 2 × period → `ArgumentException` |

## 4. Test cases (expected values = compiled TRF 4.10.0)

| ID | Input | Expected |
|---|---|---|
| T1 | TRF test_seqs s1/s2/s3, maxPeriod 2000 | README tables: 1–35 p7 5.0 ×7 100/0 70; 1–84 p12 7.0 168; 1–1225 p35 35.0 2450; composition/entropy |
| T2 | `CACACACACA`, min 10 | 1–10 p2 5.0 CA 100/0 20 |
| T3 | `CAGCAGCAGTAGCAGCAG`, p3, min 10 | 1–18 p3 6.0 CAG 86.67/0 27 (13/2/0) |
| T4 | `CACACATACACA`, min 10 | nothing (sum of heads 4 < 5) |
| T5 | 29-bp CAG deletion tract, p3 | 1–29 p3 10.0 CAG 92.59/7.41 51 (25/0/2) |
| T6 | same, maxPeriod 500, min 10 / 50 | rows p3 (51), p14 (56), p11 (44) / p3, p14 |
| T7 | flanked `TT…GACCA` + lowercase | same rows shifted (3–31, 4–31, 7–28); case-insensitive |
| T8 | `NNNNNNNN`, N×200 | nothing |
| T9 | N inside CAG array | 5–34 p3 10.0 92.59/0 51 |
| T10 | A×30 G A×26 | 1–57 p1 57.0 96.43/0 105, entropy 0.13 |
| T11 | (CA)×30, maxPeriod 500; with minPeriod 3 | one row p2 120; nothing |
| T12 | three random sequences with embedded repeats | all TRF rows (2 + 3 + 2) reproduced |
| T13 | `TrfSumOfHeadsCriterion(d)` for d = 1, 22, 23, 29, 30, 50, 100, 159, 160, 500, 2000 | TRF `sumdata80`: 5, 5, 6, 9, 6, 15, 39, 69, 43, 177, 818 |
| T14 | defaults / thresholds / validation | min 50 suppresses CA×5; deletion tract 51 reported at 50, not at 52; V-1 |
| B1–B3 | Bernoulli on CA×5, CAG/TAG tract, CA×6 with T | 8/0/0, 13/2/0 (PM 0.8667), 8/2/0 (PM 0.80, meets 0.80) |
| B4 | `ACACTGTG`, p4 | 0 trials, 0 pairs, PM 0, not meeting |
| B6 | deletion tract, p3 | 25/0/2, 27 trials, 9 pairs |
| B-link | Bernoulli on a detected region | equals the reported % matches / % indels |

### TRF parameters and outputs (WP6; `RepeatFinder_TrfParameters_Tests.cs`; crafted U1–U5, expected = compiled TRF)

| ID | Input | Expected (TRF .dat / files) |
|---|---|---|
| P1 | U1 (2 N inside a period-7 array), recommended | 61–124 p7 9.1 7 92/0 110 14/14/26/42 entropy 1.83; `EntropyTrf` = 1.829258111162015 (denominator 64); legacy `Entropy` 1.8425 (prints 1.84) |
| P2 | U3, U4 pure ACGT | `EntropyTrf = Entropy` (1.18, 1.66) |
| P3 | `Recommended` | 2 7 7 80 10 50 500, `-l` 2 000 000, redundancy on, no flanks |
| P4 | U4, recommended | 41–136 p12 8.0 12 86/4 140 TCTTCACTGCCC; 43–136 p49 score 152 |
| P5 | `2 3 5 80 10 40 200` on U3/U4/U1/U2 | U3 1–60 p2 89/0 105; U4 score 160; U1 score 118; U2 51–173 p25 5.0 24 82/6 178 and p49 2.5 50 86/2 196 |
| P6 | PM 75: `2 7 7 75 20 50 500`, `2 5 5 75 10 30 100` | U2 51–167 p25 4.8 24 84/10 148 / 168; U3 score 101; U1 score 114 |
| P7 | `TrfSumOfHeadsCriterion(d, 75)` | TRF `sumdata75`: d 1/18/20/21/29/30/43/44/100/159/160/500/2000 → 5/5/5/6/10/6/11/7/27/50/24/116/567 |
| P8 | `TrfSumOfHeadsCriterion(d, 80)` | = default overload (5, 6, 43, 818) |
| P9 | `-r` (`EliminateRedundancy = false`) | U1 periods 7/14/21 (110); U3 2/4/6 (95); U5 copies 100.0/50.0/33.3 (400) |
| P10 | `MaxRepeatLength` 150 / 60 (TRF built with `-l` in bp) | U5 79–228 75.0 300 / 169–228 30.0 120; U1 61–117 8.1 96 |
| P11 | invalid parameter sets | V-2 |
| P12 | `-m` | U1 / U3 / U4 equal to TRF's masked file; U2 soft mask = lower-cased 51..167 |
| P13 | `FlankLength` 50 / 500 | U1 flanks = TRF `-ngs` 50-bp flanks; U3 left flank empty (TRF '.'), right = rest of sequence (`-f`) |
| P14 | alignment rows, U4 | equal to TRF alignment file (`TCTTC-GTGCCC`, `TCTT-CACTGCCC` columns) |

### TRF detection pipeline (WP7; `RepeatFinder_TrfDetection_Tests.cs`; expected = compiled TRF)

Each sequence was chosen because the WP6 implementation disagreed with TRF on it and the disagreement is removed only
by the named component (ablation builds; Evidence §WP7).

| ID | Input | Expected |
|---|---|---|
| D1 | `TrfApparentSize(100, 75)` | 56 (README example) |
| D2 | `TrfApparentSizeOffset(d, pm)` at 14 (d, PM) where the exact value equals TRF's `waitdata` entry | 18, 18, 19, 28, 29, 31, 65, 67, 68 (PM 80, d 1/20/21/30/44/60/250/1000/2000); 15, 15, 25, 39, 44 (PM 75, d 1/25/30/50/150) |
| D3 | offsets for all d, both PM | INV-12 |
| D4 | apparent size, recommended set, 322-bp / 351-bp sequences | TRF reports nothing (without the test: period-28 148–209 / period-22 129–182 rows) |
| D5 | apparent size, 723-bp period-267 array | 65–658 p267 2.2 264 85/11 788 (without: consensus 265, score 790) |
| D6 | narrow band, period 31 / 36 / 264 arrays (266 / 279 / 760 bp) | 102–182 p31 2.6 30 84/3 99; 67–165 p36 2.8 37 92/6 168; 96–691 p264 2.3 263 83/9 761 (full WDP: 102–162 / 95, consensus 36 / 164, 683 / 752) |
| D7 | narrow band via the legacy overload (maxPeriod 500) | same row as D6 period 31 |
| D8 | active distances, `2 5 5 75 10 30 100`, 175 bp | single row 1–175 p43 4.0 42 83/9 256 (all distances active: extra period-23 row) |
| D9 | best-period list (d > 250), `2 7 7 80 10 50 2000`, 1 374 bp | 58–1201 p278 4.1 275 79/10 1140 (without the list: 36–1201, score 1185) |

## 5. Cross-check / Differential Oracle

- Per-candidate analysis vs instrumented TRF: 1 524 / 1 524 identical (pattern ≤ 20), 827 / 959 (> 20, TRF band).
- Whole pipeline vs TRF `.dat` (700 random sequences): 92.6 % exact rows / 96.0 % region level (periods ≤ 20);
  80.5 % / 93.1 % (periods ≤ 100). Details in the Evidence file.
- WP6 (new 700-sequence set, 25 % with N; `TandemRepeatsFinderParameters` API): recommended 2 7 7 80 10 50 500 →
  87.8 % exact / 99.4 % region (1 305 TRF rows; periods ≤ 20: 98.3 % / 99.8 %); six non-default sets 84–92 % /
  99.1–99.9 % except the most permissive 2 3 3 80 20 (64.2 % / 93.8 %); `-r`: 89.1 % / 99.7 %; `EntropyTrf` equal to
  TRF on every same-locus N-containing row (229/229 … 566/566); masks identical whenever the loci are identical
  (601/601, 565/565, 582/582, 461/461 sequences); flanks 50 bp 100 %, 500 bp 314/314; alignment rows 191/191
  (consensus ≤ 20), 116/123 (> 20).

- WP7 (TRF detection pipeline complete; same 700 sequences): exact / region = 2 7 7 80 10 50 500 99.8 / 100
  (1 303/1 305); 2 5 7 80 10 50 2000 99.9 / 100; 2 3 5 80 10 40 200 99.9 / 100; 2 7 7 75 20 50 500 99.9 / 100;
  2 5 5 75 10 30 100 99.8 / 100; 3 7 7 80 10 60 50 100 / 100; 2 3 3 80 20 50 500 99.8 / 100; legacy overload
  (maxPeriod 500) 99.8 / 100; `-r` (two sets) 99.8 % / 99.9 % exact, 100 % region; `-l` 60 / 120 / 250 bp 100 / 100 / 99.8 %; masks 700/700,
  700/700, 700/700, 699/700 identical; alignment rows 192/192 + 149/149; 1 Mb sequence 1 543/1 544. With TRF's own
  simulated apparent-size table substituted (diagnostic only) every one of these is 100 %.

## 6. Declared residual

None in the method. The apparent-size cut-offs are the exact values of the distribution TRF estimates by Monte-Carlo
simulation; TRF's shipped tables carry simulation noise (equal at 825/2000 and 713/2000 distances, |Δ| ≤ 3 / 4).
Each of the 25 remaining non-identical (parameter set, sequence) cases (16 TRF rows) over the seven sets traces (bisection) to a single table
entry where TRF's noisy value is 1 below the exact one (Evidence §WP7). TRF exits when a band exceeds 150 cells (PI 20,
patterns ≥ 1365); this implementation has no such limit. `Entropy` stays the normalised Shannon entropy; `EntropyTrf`
is TRF's column in every case (WP6).
