# Test Specification: REP-APPROX-001

**Test Unit ID:** REP-APPROX-001
**Area:** Repeats
**Algorithm:** Approximate (TRF) Tandem-Repeat Detection + TRF Bernoulli statistics
**Status:** ☑ Complete — re-validated 2026-09-30 (campaign 2026-09, batch B04; Stage A corrected, Stage B fixed)
**Last Updated:** 2026-09-30

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

## 2. Canonical Method(s)

| Method | Type |
|---|---|
| `FindApproximateTandemRepeats(DnaSequence, minPeriod = 1, maxPeriod = 6, minScore = 50)` | Canonical |
| `FindApproximateTandemRepeats(string, …)` | Overload (case-insensitive; null/empty → empty) |
| `ComputeBernoulliStatistics(string repeatTract, int period, double expectedMatchProbability = 0.80)` | Canonical |
| `TrfSumOfHeadsCriterion(int d)` (internal) | Helper (tested) |

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_ApproximateTandemRepeats_Tests.cs`
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

## 5. Cross-check / Differential Oracle

- Per-candidate analysis vs instrumented TRF: 1 524 / 1 524 identical (pattern ≤ 20), 827 / 959 (> 20, TRF band).
- Whole pipeline vs TRF `.dat` (700 random sequences): 92.6 % exact rows / 96.0 % region level (periods ≤ 20);
  80.5 % / 93.1 % (periods ≤ 100). Details in the Evidence file.

## 6. Declared residual

TRF's apparent-size criterion (simulated), random-walk range distances, best-period list (d > 250) and narrow-band
WDP (patterns > 20) are not reproduced; entropy is the canonical normalised Shannon entropy (differs from TRF only
for regions containing N).
