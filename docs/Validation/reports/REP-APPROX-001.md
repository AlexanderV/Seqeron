# Validation Report: REP-APPROX-001 — Approximate (TRF) Tandem-Repeat Detection

> **Re-validation 2026-09 / 2026-10 (campaign 2026-09, batch B04) — supersedes the 2026-06 verdict below.**
> The 2026-06 review had no TRF binary; TRF 4.10.0 compiled from source showed the old statistics and detection were
> not TRF's (B04 F14–F18: analysis now TRF-identical for patterns ≤ 20, k-tuple + sum-of-heads detection).
> Completeness audit WP6 (B04 F40–F42): `TandemRepeatsFinderParameters` overloads expose every TRF parameter
> (weights, PM 80/75 with both sum-of-heads tables reproduced 2000/2000, PI random-walk range, Minscore, MaxPeriod /
> MAXDISTANCE, `-l`, `-r`, `-f`); `MaskApproximateTandemRepeats` (`-m`, + soft mask); `EntropyTrf`, alignment rows and
> flanks on the result. Cross-check vs compiled TRF on 700 sequences: recommended set 87.8 % rows exact / 99.4 %
> region level, non-default sets comparable (Evidence, WP6 revision); masks 0 mismatches given identical loci;
> `EntropyTrf` 0 mismatches on every N-containing row.
> Completeness audit WP7 (B04 F43–F46): the TRF detection pipeline is complete — apparent-size criterion (the
> distribution TRF simulates, derived exactly; README example y = 56 reproduced; equal to TRF's simulated table at
> 825/2000 and 713/2000 distances, within its noise elsewhere), random-walk range sums over active distances only,
> best-period list for d > 250, narrow-band WDP for patterns > 20 (also on the legacy overloads, which now take the
> range too). Compiled TRF on the same 700 sequences: **99.8–100 % rows exact, 100 % region level on all seven
> parameter sets** (recommended 1 303/1 305; permissive 2 3 3 80 20 64.2 % → 99.8 %), `-r` 99.8/99.9 %, `-l` 100 /
> 100 / 99.8 %, masks 699–700/700, alignment rows 149/149 for consensus > 20, 1 Mb sequence 1 543/1 544; 3× faster.
> Each of the 25 remaining (parameter set, sequence) cases (16 TRF rows) is traced (bisection) to one entry of TRF's Monte-Carlo apparent-size table that is
> 1 below the exact value; with TRF's table substituted every comparison is 100 %. **Stage A: corrected · Stage B:
> fixed · State: FIXED** — no method residual (TRF's simulation noise and its 150-cell band abort are documented).
> Current sources: `tests/TestSpecs/REP-APPROX-001.md`, `docs/Evidence/REP-APPROX-001-Evidence.md`,
> `docs/Validation/review-2026-09/B04.md` (F14–F18, F40–F46). Tests:
> `Unit/Analysis/RepeatFinder_ApproximateTandemRepeats_Tests.cs`, `Unit/Analysis/RepeatFinder_TrfParameters_Tests.cs`,
> `Unit/Analysis/RepeatFinder_TrfDetection_Tests.cs`.

---

*Historical 2026-06 report (superseded):*

- **Validated:** 2026-06-25   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindApproximateTandemRepeats` (string + `DnaSequence` overloads), `RepeatFinder.ComputeBernoulliStatistics`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
`FindApproximateTandemRepeats`, `ComputeBernoulliStatistics`

- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs` (lines ~280–671, +result records ~1059–1090)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_ApproximateTandemRepeats_Tests.cs`

## Authoritative sources opened (this session)
1. **Benson G (1999)** "Tandem repeats finder: a program to analyze DNA sequences", *Nucleic Acids Res* 27(2):573–580, https://doi.org/10.1093/nar/27.2.573 (PMC148217). Confirms: alignment weights "**2 for match, 7 for mismatch and indel** weights as suggested values"; the Bernoulli model "**we model the alignment of two tandem copies of a pattern of length n by a sequence of n independent Bernoulli trials**"; consensus determined "**by majority rule**"; sum-of-heads statistic R and a Minscore threshold for reporting.
2. **TRF definitions page** (tandem.bu.edu/trf/trf.definitions.html). Confirms verbatim: Period Size = estimated pattern size; Copy Number = "Number of copies aligned with the consensus pattern"; Consensus Size = "Size of consensus pattern (may differ slightly from the period size)"; **Percent Matches = "Percent of matches between adjacent copies overall"**; **Percent Indels = "Percent of indels between adjacent copies overall"**; Alignment Score uses Smith–Waterman style with match +2, mismatch/indel ∈ {3,5,7} negative; and critically **"Statistics refers to the matches, mismatches and indels overall between adjacent copies in the sequence, not between the sequence and the consensus pattern."**
3. **TRF parameters / README** — defaults **PM = .80 and PI = .10** ("PM=80 and PI=10"); probabilistic data available for PM∈{80,75}, PI∈{10,20}, best with PM=80/PI=10.

## Stage A — Description

| Checklist item | Result |
|---|---|
| Source quality | ✅ Primary literature (Benson 1999) + official TRF docs retrieved this session; claims read directly, not from citation labels. |
| Formula correctness | ✅ Weights +2/−7/−7 are the Benson "2 7 7" suggested set. Bernoulli model, majority-rule consensus, Minscore reporting threshold all match the source. |
| Definitions & conventions | ✅ Copy number, consensus size, percent matches, percent indels match the TRF definitions page verbatim; %matches/%indels are over alignment columns. |
| Edge-case semantics | ✅ "two or more contiguous copies" ⇒ ≥2 copies required (single char / period-too-long ⇒ none; Bernoulli needs ≥2 copies). N handled as a literal base (defined; no IUPAC special-casing). |
| Independent cross-check | ✅ Hand-computed on controlled tracts (below). TRF binary not installable (`which trf` = none, no conda); hand-computation + Benson's worked definitions are the oracle, as the protocol permits. |
| Invariants | ✅ %matches ∈ [0,100]; reported score ≥ Minscore; match+mismatch+indel partition the columns; deterministic. |

### Independent hand-computation (oracle = Benson 1999), reproduced by the C#

| Case | Input | Consensus | Match / Mis / Indel | Score (+2/−7/−7) | %matches | %indels | Copies |
|---|---|---|---|---|---|---|---|
| A1 perfect | `CACACACACA` p=2 | `CA` | 10 / 0 / 0 | **20** | **100** | 0 | 5 |
| A2 1 sub | `CAGCAGCAGTAGCAGCAG` p=3 | `CAG` | 17 / 1 / 0 | **27** | **94.444…** | 0 | 6 |
| A3 1 sub | `CACACATACACA` p=2 | `CA` | 11 / 1 / 0 | **15** | **91.667…** | 0 | 6 |
| A4 1 del | `CAGCAGCAGCAGCAGAGCAGCAGCAGCAG` (29 bp) p=3 | `CAG` | 29 / 0 / 1 | **51** | **96.667…** | **3.333…** | 29/3 |

Bernoulli (adjacent-copy) cross-check:

| Case | Copies | Pairs | Trials | Match/Mis/Indel | PM | PI |
|---|---|---|---|---|---|---|
| B1 perfect CA×5 | CA·CA·CA·CA·CA | 4 | 8 | 8/0/0 | **1.0** | 0 |
| B2 CAG, copy4=TAG | …·TAG·… | 5 | 15 | 13/2/0 | **13/15** | 0 |
| B3 CA, idx6=T | …·TA·… | 5 | 10 | 8/2/0 | **0.80** | 0 |
| B4 ACAC/TGTG | 2 | 1 | 4 | 0/4/0 | **0.0** | 0 |

All values above were produced by running the fixture and confirmed equal to the hand-computation.

**Stage A verdict: ✅ PASS.** Description (weights, statistics definitions, Bernoulli model, Minscore, majority-rule consensus) is faithful to Benson 1999 and the TRF docs.

## Stage B — Implementation

- **Code path reviewed:** `RepeatFinder.cs` — `TrfScoring` (Match +2, Mismatch −7, GapOpen/GapExtend −7); `FindApproximateTandemRepeatsCore` (per-start/per-period window growth, majority consensus, tile-to-whole-copies reference, global align, non-overlap suppression by score); `MajorityConsensus`; `TileTo`; `ComputeTrfStatistics` (column counts, %matches/%indels over columns, copy number = non-gap bases / period); `ComputeBernoulliStatistics` (adjacent-copy pairwise alignment, match/mismatch/indel partition, PM/PI, E[matches]=PM·d, ≥-threshold flag).
- **Flat-indel realisation:** verified `SequenceAligner.GlobalAlignCore` uses a **linear** gap (`GapExtend` only; `GapOpen` is unused in the global DP), so the −7 per gap column matches TRF's flat indel weight exactly. (Confirmed at `SequenceAligner.cs:236–242, 256–257`.)
- **Cross-verification:** all A1–A4 and B1–B4 numbers above reproduced by the code (fixture green).
- **Variant/delegate consistency:** `DnaSequence` overload delegates to the string core (test `…DnaSequenceOverload_MatchesStringOverload`); determinism asserted (`…RunTwice_ProducesIdenticalResults`).
- **Numerical robustness:** integer column counts; %·100 over column count guarded for 0 columns; no overflow on stated ranges.

### Test-quality audit (HARD gate)
- Expected values trace to Benson 1999 / hand-computation, **not** code echoes (literal `20`, `27`, `15`, `51`, `13/15`, `0.80`, etc., each derived from the alignment columns in comments).
- Every public method/overload covered: both `FindApproximateTandemRepeats` overloads, `ComputeBernoulliStatistics` (default + custom PM), the exposed constants (`DefaultApproximateMinScore`, `TrfDefaultMatchProbability`, `TrfDefaultIndelProbability`).
- Stage-A paths/edge cases covered: perfect (100%/0%), substitution (lowers %matches), indel path (%indels>0), Minscore threshold (suppress 20<50, report 51≥50), ≥2-copy requirement (Bernoulli throws on 1 copy), empty, no-repeat, invalid period, null, single-char, period-longer-than-seq, all-N (literal). "minReps 0" is N/A: this API hard-codes the ≥2-copy minimum (no minReps parameter) — the analogous boundary is covered.
- No green-washing: assertions check exact sourced values, exceptions for invalid inputs, deterministic.

### Changes made this session
- Added 3 edge-case tests to `RepeatFinder_ApproximateTandemRepeats_Tests.cs` (A12 single-char ⇒ empty; A13 period-longer-than-seq ⇒ empty; A14 all-N ⇒ literal perfect homopolymer, score 16/100%/0%/8 copies). No production code changed.

**Stage B verdict: ✅ PASS.** Code faithfully realises the validated description; no defect.

## Documented scope boundary (not a defect)
Genome-scale **k-tuple seeding** (the R(d,k,PM) sum-of-heads percentile cut-offs and the W(d,PI) random-walk distance band, which rely on TRF's non-redistributable simulation tables) is **not** reproduced; the detector uses deterministic exhaustive per-start/per-period window evaluation instead. This is a documented performance boundary, not a correctness gap — the per-repeat statistics (consensus, %matches, %indels, score, PM/PI) are exact.

## Verdict & follow-ups
- **Stage A: ✅ PASS · Stage B: ✅ PASS · State: ✅ CLEAN.** No defect found; coverage strengthened with 3 sourced edge-case tests. Full unfiltered `dotnet test Seqeron.sln -c Debug` green (Seqeron.Genomics.Tests 18778 passed, 0 failed). No follow-ups.
