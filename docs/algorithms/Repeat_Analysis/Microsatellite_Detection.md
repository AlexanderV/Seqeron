# Microsatellite Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-STR-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-09-30 |

## 1. Overview

Microsatellite detection identifies short tandem repeats (STRs, also called microsatellites or simple sequence repeats) whose motif length is between 1 and 6 nucleotides [1][3][4]. The repository implements exact consecutive-repeat detection in `RepeatFinder.FindMicrosatellites` (the default), classifies each hit by repeat-unit length, and exposes overloads for `DnaSequence`, raw strings, and cancellation-aware execution. The implementation also removes redundant compound motifs such as `ATAT` when they are just repetitions of a smaller motif, and it suppresses results fully contained inside already reported repeat intervals. An **opt-in approximate detector**, `RepeatFinder.FindApproximateTandemRepeats`, additionally finds **imperfect / interrupted** tandem repeats (those containing substitutions or indels) using the Tandem Repeats Finder (TRF) model [6][7]: k-tuple matches at a common distance that pass Benson's sum-of-heads criterion trigger a wraparound-dynamic-programming (WDP) alignment of the sequence against tandem copies of the candidate pattern, the consensus is taken by majority rule, the sequence is realigned against it, and the TRF table (indices, period, copy number, consensus size, % matches and % indels between adjacent copies, score, composition, entropy) is reported (review 2026-09, REP-APPROX-001: identical to compiled TRF 4.10.0 on every analysed candidate with pattern ≤ 20). The default perfect-repeat detector is unchanged. Microsatellite biology matters clinically, forensically, and evolutionarily because repeat expansions drive many genetic disorders and STR polymorphism underpins DNA profiling [1][2][3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Microsatellites are tandemly repeated DNA motifs 1-6 bp long, typically repeated about 5-50 times in natural genomic settings [1][3][4]. Standard terminology is:

| Class | Unit Length | Example |
|-------|-------------|---------|
| Mononucleotide | 1 bp | `AAAAAA` |
| Dinucleotide | 2 bp | `CACACACA` |
| Trinucleotide | 3 bp | `CAGCAGCAG` |
| Tetranucleotide | 4 bp | `GATAGATAGATAGATA` |
| Pentanucleotide | 5 bp | any repeated 5-bp motif |
| Hexanucleotide | 6 bp | any repeated 6-bp motif |

Trinucleotide repeat expansions are associated with more than 30 genetic disorders [2]. The legacy reference set highlights Huntington's disease, Fragile X syndrome, Friedreich's ataxia, and myotonic dystrophy as canonical examples [2]. Forensic DNA profiling also depends on STR variability, with tetra- and pentanucleotide repeats commonly preferred because they reduce PCR stutter relative to shorter motifs [1].

### 2.2 Core Model

For a motif $U$ of length $m$ and repeat count $k$, a microsatellite occupies a contiguous region when:

$$
S[p..p + km) = U^k \quad \text{with} \quad 1 \le m \le 6
$$

The implementation searches candidate motif lengths from `minUnitLength` through `maxUnitLength`, skips motifs that are themselves repetitions of a smaller motif, counts consecutive copies, and emits a result when the count reaches `minRepeats`. Each result is classified into a `RepeatType` by motif length.

**Approximate (TRF) model.** Benson (1999) defines a tandem repeat as "two or more contiguous, *approximate* copies of a pattern of nucleotides" and reports, for each repeat, the period size, the number of copies aligned with the consensus pattern, the consensus size, the percent of matches and percent of indels between adjacent copies overall, and an alignment score [6]. The alignment is scored Smith-Waterman style with weights for match, mismatch and indels; the recommended parameter set is match `+2`, mismatch `7`, indel (delta) `7` (applied as negatives), and only repeats scoring at least `Minscore = 50` are reported [6]. The consensus pattern is determined "by majority rule from the alignment" [6]. The alignment is a *local* wraparound DP: the sequence (rows) is aligned against unlimited tandem copies of the pattern (columns wrap from the last pattern position to the first), each row computed in two passes; the score is the best local alignment score [6][7]. The statistics compare each copy with the next one *through* the consensus alignment ("between adjacent copies in the sequence, not between the sequence and the consensus pattern" [7]); the period is "the most common matching distance between corresponding characters in the alignment" and may differ from the consensus size [7].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every result satisfies `RepeatCount >= minRepeats`. | Results are emitted only after counting the required number of consecutive copies. |
| INV-02 | `minUnitLength <= RepeatUnit.Length <= maxUnitLength`. | The outer search loop enumerates only those unit lengths. |
| INV-03 | `TotalLength = RepeatUnit.Length × RepeatCount`. | `MicrosatelliteResult.TotalLength` is constructed from those two fields. |
| INV-04 | `RepeatType` matches the reported unit length. | `ClassifyRepeatType` maps unit lengths 1 through 6 to the corresponding repeat class. |
| INV-05 | Per unit length, each maximal perfect run is reported exactly once, at its left end: `Position = 0` or `S[Position−1] ≠ S[Position−1+p]` (left-maximal), and `RepeatCount = ⌊runLength/p⌋` (right-maximal; a trailing partial copy is not counted). Rotations of the same run (e.g. `TA` inside `ATATATA`) are never reported. `RepeatUnit` is primitive and consists of A/C/G/T only. | The scan visits only left-maximal positions, extends the run while `S[x] = S[x−p]`, and skips to `e−p+1` (Kolpakov & Kucherov 1999 maximal repetitions [8]; MISA leftmost match [9]; pytrf/Krait run start [10]). |
| INV-06 | Every approximate result has `AlignmentScore >= minScore` and `CopyNumber >= 1.9` (`>= 1.8` for consensus > 100). | TRF report rules [6][7]. |
| INV-07 | For an approximate result, `PercentMatches + PercentIndels <= 100` (both in [0,100]); a perfect tract yields 100 / 0. | Matches, mismatches and indels partition the adjacent-copy comparisons [7]. |
| INV-08 | For an approximate result, `CopyNumber` = aligned consensus columns / `ConsensusSize`; `Period` = most common distance between matching characters of adjacent copies (may differ from `ConsensusSize`); `|Consensus| = ConsensusSize`. | TRF "Table Explanation" and "Consensus Pattern and Period Size" [7]. |
| INV-08b | No two reported approximate repeats overlapping by ≥ 90 % of one of them are redundant (same period and no higher score, or a multiple period scoring ≤ 1.1×). Output ordered by `Start`. | TRF redundancy elimination [7]. |

| INV-09 | (Bernoulli) For `ComputeBernoulliStatistics`, `MatchProbability ∈ [0,1]`, `IndelProbability ∈ [0,1]`, and `Matches + Mismatches + Indels = BernoulliTrials`; a perfect tract yields `MatchProbability = 1`, `IndelProbability = 0`; the values equal TRF's adjacent-copy counts for the tract. | Each Bernoulli trial is one adjacent-copy comparison of the TRF consensus alignment, classified as exactly match / mismatch / indel [6][7]. |
| INV-10 | (Bernoulli) `ExpectedMatches = MatchProbability × BernoulliTrials` and `MeetsExpectedMatchProbability ⇔ MatchProbability ≥ expectedMatchProbability`. | `ExpectedMatches` is the Bernoulli mean E[heads] = PM·d; the flag is the direct comparison to the assessed PM [6]. |

> A = perfect STR detection (`FindMicrosatellites`), B = approximate / imperfect tandem-repeat detection (`FindApproximateTandemRepeats`, TRF model [6]), C = TRF Bernoulli statistical measures (`ComputeBernoulliStatistics`, Benson 1999 [6]). Invariants INV-01..INV-05 govern A; INV-06..INV-08 govern B; INV-09..INV-10 govern C.

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to search. | `DnaSequence` overloads throw on `null`; string overloads yield no results for `null` or empty input. |
| `minUnitLength` | `int` | `1` | Minimum repeat-unit length. | Values below `1` throw `ArgumentOutOfRangeException`. |
| `maxUnitLength` | `int` | `6` | Maximum repeat-unit length. | Values below `minUnitLength` throw `ArgumentOutOfRangeException`. |
| `minRepeats` | `int` | `3` | Minimum number of consecutive copies to report. | Values below `2` throw `ArgumentOutOfRangeException`. |
| `cancellationToken` | `CancellationToken` | optional | Cancellation support for long-running scans. | Used only by the cancellable overloads. |
| `progress` | `IProgress<double>?` | optional | Progress callback for cancellable scans. | Receives values from `0.0` to `1.0` in the cancellable implementation. |
| `minPeriod` | `int` | `1` | (approximate) Minimum reported period; applied after redundancy elimination (never resurrects a redundant multiple). | `FindApproximateTandemRepeats`; values below `1` throw `ArgumentOutOfRangeException` (eager, both overloads). |
| `maxPeriod` | `int` | `6` | (approximate) Maximum candidate distance and reported period (TRF `MaxPeriod`; TRF recommends 500). | Values below `minPeriod` or above `2000` (TRF 4.10.0 limit) throw `ArgumentOutOfRangeException`. |
| `minScore` | `int` | `50` | (approximate) Minimum TRF alignment score to report. | Values below `1` throw `ArgumentOutOfRangeException`; default `DefaultApproximateMinScore = 50` per Benson (1999) [6]. |
| `repeatTract` | `string` | required | (Bernoulli) Observed tandem-repeat tract (≥ 2 copies). | `ComputeBernoulliStatistics`; `null` throws `ArgumentNullException`; fewer than two copies throws `ArgumentException`. |
| `period` | `int` | required | (Bernoulli) Candidate period; the last `period` bases are the initial pattern. | `ComputeBernoulliStatistics`; values outside `1..2000` throw `ArgumentOutOfRangeException`. |
| `expectedMatchProbability` | `double` | `0.80` | (Bernoulli) PM threshold the tract is assessed against. | `ComputeBernoulliStatistics`; values outside `[0,1]` (or NaN) throw `ArgumentOutOfRangeException`; default `TrfDefaultMatchProbability = 0.80` per Benson (1999) [6]. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | 0-based start position of the reported microsatellite. |
| `RepeatUnit` | `string` | Repeated motif. |
| `RepeatCount` | `int` | Number of consecutive copies of the motif. |
| `TotalLength` | `int` | Total length of the repeat tract in bases. |
| `RepeatType` | `RepeatType` | Unit-length classification from mono- through hexanucleotide. |

`FindApproximateTandemRepeats` returns `ApproximateTandemRepeatResult`:

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | 0-based start of the repeat (TRF prints 1-based: TRF start = `Start + 1`, TRF end = `Start + SpanLength`). |
| `SpanLength` | `int` | Length of the repeat region (all sequence symbols between its first and last aligned base). |
| `Period` | `int` | TRF period: most common distance between matching characters of adjacent copies [7]. |
| `ConsensusSize` | `int` | Size of the consensus pattern (may differ from `Period`) [7]. |
| `Consensus` | `string` | Majority-rule consensus, starting at the phase of the first repeat base [6][7]. |
| `CopyNumber` | `double` | Copies aligned with the consensus = aligned consensus columns / `ConsensusSize` [7]. |
| `PercentMatches` | `double` | Percent of matches between adjacent copies overall (exact; TRF truncates to an integer) [7]. |
| `PercentIndels` | `double` | Percent of indels between adjacent copies overall (exact) [7]. |
| `AlignmentScore` | `int` | WDP local alignment score against the consensus [6]. |
| `PercentA` / `PercentC` / `PercentG` / `PercentT` | `double` | Composition of the region (denominator = `SpanLength`, so N lowers all four) [7]. |
| `Entropy` | `double` | Shannon entropy (bits, 0–2) of the region's A/C/G/T composition, via canonical `SequenceComplexity.CalculateShannonEntropy` [7]. |

`ComputeBernoulliStatistics` returns `TandemRepeatBernoulliStatistics` (the TRF Bernoulli statistical measures, Benson 1999 [6]):

| Field | Type | Description |
|-------|------|-------------|
| `Period` | `int` | Repeat period (copy length). |
| `AdjacentCopyPairs` | `int` | ⌈aligned copy number⌉ − 1 (0 when fewer than two copies align). |
| `BernoulliTrials` | `int` | Adjacent-copy comparisons of the TRF consensus alignment [6][7]. |
| `Matches` / `Mismatches` / `Indels` | `int` | Outcome counts (heads = matches) — equal to TRF's "Matches / Mismatches / Indels" statistics for the tract. |
| `MatchProbability` | `double` | PM = P(Heads) = average percent identity between adjacent copies, as a fraction (0–1) [6]. |
| `IndelProbability` | `double` | PI = average percentage of insertions/deletions between adjacent copies, as a fraction (0–1) [6]. |
| `PercentMatches` / `PercentIndels` | `double` | PM / PI as percentages (0–100) [6]. |
| `ExpectedMatches` | `double` | Bernoulli mean E[heads] = `MatchProbability × BernoulliTrials` [6]. |
| `MeetsExpectedMatchProbability` | `bool` | `true` when `MatchProbability ≥ expectedMatchProbability` (default PM = 0.80, Benson 1999 [6]). |

### 3.3 Preconditions and Validation

All overloads validate `minUnitLength`, `maxUnitLength`, and `minRepeats`, rejecting `minUnitLength < 1`, `maxUnitLength < minUnitLength`, and `minRepeats < 2`. `DnaSequence` overloads throw `ArgumentNullException` for `null` sequence input. String overloads yield no results for `null` or empty input and normalize non-empty input to uppercase before scanning. Reported positions are 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize raw-string input to uppercase when needed.
2. For each unit length `p` from `minUnitLength` to `maxUnitLength`, scan positions `i` where at least `minRepeats` copies could fit.
3. Skip `i` unless it is a run start (left-maximal): `i = 0` or `S[i−1] ≠ S[i−1+p]`.
4. Extend the run to its right-maximal end `e` (exclusive) while `S[e] = S[e−p]`; `RepeatCount = ⌊(e−i)/p⌋`.
5. If `RepeatCount ≥ minRepeats` and the unit `S[i..i+p)` is primitive (not a power of a shorter word — MISA "reject false type motifs") and contains only A/C/G/T (MISA `[acgt]`, pytrf skips `N`), emit one `MicrosatelliteResult`.
6. Continue at `e−p+1` (every position in `(i, e−p]` lies inside the same run). O(n) comparisons per unit length.

For the approximate detector (`FindApproximateTandemRepeats`, TRF model [6][7]):

1. Scan positions `i` (1-based internally, as TRF) and distances `d = 1..maxPeriod`. A *k-tuple match* at distance `d` ends at `i` when the last `k` symbols equal those `d` earlier (A/C/G/T only); `k = 4 / 5 / 7` for `d ≤ 29 / 30..159 / ≥ 160` (Benson 1999 Table 1, PM = .80). Runs of adjacent tuple matches are kept in a window of the last `max(d, 20)` positions.
2. **Sum-of-heads criterion:** the heads in the window must reach `max(k+1, ⌊μ − 1.65σ⌋)`, `μ, σ` = exact mean / s.d. of `R(d,k,PM)` (heads in runs ≥ k of a Bernoulli(0.80) sequence of length `d`) — reproduces TRF's `sumdata80` table for all `d ≤ 2000`. Positions already covered by an alignment at the same `d` are skipped.
3. **WDP** with pattern `S[i−d+1..i]`: a backward local scan finds the leftmost row reaching the best score; a forward local alignment starting one pattern length earlier gives the optimum (zero cells beyond `i` / before `i − max(d,20)` are killed so the alignment stays contiguous with the candidate). Traceback preference: match/mismatch, then sequence-vs-gap, then pattern-vs-gap.
4. Require ≥ 1.9 copies (ramp to 1.8 for 50–100, 1.8 above 100) and `d` among the **three best periods** of the aligned region (dinucleotide-distance histogram with its least-squares trend removed; period 1: ≥ 80 % one base).
5. **Consensus** by majority rule (a position is deleted when gaps are at least as frequent; an insertion point receives its most frequent base when insertions occur in ≥ half the passes), then **realign** against the consensus (step 3) and re-apply the copy rule and `minScore`.
6. **Statistics** from the final alignment: two cursors one consensus period apart compare each copy with the next (match / mismatch / indel), period = most common matching distance, composition and entropy over the region.
7. Drop periods above `maxPeriod`, sort by start, apply TRF **redundancy elimination** (≥ 90 % overlap: same period with no higher score, or a multiple period with score ≤ 1.1×, is removed), drop periods below `minPeriod`, order by (start, end, period).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Approximate-detector scoring constants (Tandem Repeats Finder recommended set "2 7 7 … 50 …" [6]):

| Constant | Value | Source |
|----------|-------|--------|
| Match weight | `+2` | Benson (1999): match weight "+2 in all options" [6] |
| Mismatch penalty | `−7` | Benson (1999): recommended Mismatch = 7 [6] |
| Indel (delta) penalty | `−7` per gap column | Benson (1999): recommended Delta = 7; flat per-column indel [6] |
| Non-ACGT symbols | never match (score −7 against anything) | TRF similarity matrix: only identical A/C/G/T score +2 ("avoid N matching itself") [7] |
| `DefaultApproximateMinScore` | `50` | Benson (1999): "Only those repeats scoring at least 50 … are reported" [6] |
| Tuple sizes | `4` (d ≤ 29), `5` (30–159), `7` (≥ 160) | Benson (1999) Table 1 / TRF 4.10.0 for PM = 80 [6][7] |
| Min distance window | `20` | TRF `Min_Distance_Window` [7] |
| Minimum copies | `1.9` (≤ 50), `1.9 − 0.002(d−50)` (50–100), `1.8` (> 100) | TRF 4.10.0 [7] |
| `MaxApproximatePeriod` | `2000` | TRF README: MaxPeriod above 2000 is an error [7] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Microsatellite detection (perfect) | `O(n × U × R)` | `O(k)` | `U` is the searched unit-length range, `R` is the average repeat count encountered while extending motifs, and `k` is the number of retained intervals/results. |
| Approximate detection (TRF) | `O(n × P + Σ_candidates (L + d) × d)` | `O(P² + L × c)` | `P = maxPeriod`; each candidate costs one WDP over its region `L` (extent-only pass, no matrix); only candidates passing the copy / best-period tests store an `L × c` traceback matrix (`c` = consensus size). Measured (Release, this container): 100 kb random, `maxPeriod` 500 → 0.34 s; 100 kb of mixed embedded repeats → 1.3 s (545 rows; TRF 549 rows); pathological 100 kb perfect `(CA)n` with `maxPeriod` 500 → 54 s (compiled TRF 74 s — every even `d` aligns the whole array in both). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable overload with progress reporting.
- `RepeatFinder.FindMicrosatellites(string, int, int, int)`: Raw-string overload with uppercase normalization.
- `RepeatFinder.FindMicrosatellites(string, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable raw-string overload.
- `RepeatFinder.FindApproximateTandemRepeats(DnaSequence | string, int minPeriod, int maxPeriod, int minScore)`: Opt-in approximate / imperfect / interrupted tandem-repeat detector (TRF model [6][7]); wraparound DP is its own algorithm (the library's pairwise aligners have no wraparound mode); entropy via canonical `SequenceComplexity.CalculateShannonEntropy`.
- `RepeatFinder.ComputeBernoulliStatistics(string repeatTract, int period, double expectedMatchProbability = 0.80)`: TRF Bernoulli measures (Benson 1999 [6]) — PM, PI and `ExpectedMatches = PM × trials` (= matches), computed **between adjacent copies** through the same TRF WDP + consensus alignment of the tract; flags whether PM ≥ the assessed PM (default 0.80).
- `RepeatFinder.TrfSumOfHeadsCriterion(int d)` (internal): Benson's sum-of-heads cut-off from the exact moments of `R(d,k,PM)`.

### 5.2 Current Behavior

Per unit length, each maximal perfect run is reported once at its left end with the run-start motif phase and the number of complete copies (review 2026-09, REP-STR-001). Redundant (non-primitive) units such as `ATAT` or `CAGCAG` are not reported — their runs appear at the primitive unit length. Units containing non-ACGT symbols (e.g. `N` in the raw-string overload) are never reported. Runs of different unit lengths may overlap and are reported independently (MISA per-motif-size convention); two distinct runs of the same period may overlap by fewer than `p` bases (e.g. `ACACACGCGCGC` → `AC×3@0`, `CG×3@5`), where greedy tools (MISA regex `/g`, pytrf `next_start`) instead restart after the first run's full copies (`GC×3@6`).

**Cross-check (2026-09):** 3,132 random + crafted cases (unit lengths 1–6, `minRepeats` 2–5, alphabets incl. `N`) agree 3,132/3,132 with an independent brute-force maximal-repetition reference; MISA-regex agreement 3,099/3,132 and pytrf 1.5.0 `STRFinder` (single motif size) 2,752/3,132, every disagreement being one of the two documented conventions (same-period runs overlapping by < p: 17; greedy tools consuming a non-primitive/`N` region first: MISA 16, pytrf 363 — pytrf has no primitivity check).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of 1-6 bp short tandem repeats (STRs) [1][3][4].
- Classification of repeats into mono-, di-, tri-, tetra-, penta-, and hexanucleotide categories [1][3].
- Consecutive-copy counting with explicit reporting of repeat unit, count, start position, and total length.

- (Approximate, TRF [6][7]) Analysis component: WDP local alignment, majority consensus, realignment, TRF statistics (indices, period, copies, consensus size, adjacent-copy % matches / % indels, score, composition), minimum copy number, three-best-periods test, redundancy elimination. **Cross-check (2026-09):** on 2 483 TRF-detected candidates (instrumented TRF 4.10.0, same position and distance) the C# analysis is identical (indices, score, consensus size, copy number, adjacent-copy counts) for all 1 524 candidates with pattern ≤ 20 (TRF's full-WDP range) and for 827 / 959 larger ones (TRF aligns those in a narrow band).
- (Approximate, TRF [6]) Detection: k-tuple trigger with Benson's tuple sizes and the sum-of-heads criterion (exact moments + normal approximation, = TRF table for all d ≤ 2000).
- (Approximate, TRF [6]) Recommended scoring constants match `+2`, mismatch `−7`, indel `−7`, and the `Minscore = 50` report threshold.

- (Bernoulli, TRF [6]) The probabilistic measures: "We model alignment of two tandem copies … by … independent Bernoulli trials"; PM = P(Heads) = "the average percent identity between the copies"; PI = "the average percentage of insertions and deletions between the copies"; statistics "between adjacent copies … not between the sequence and the consensus pattern"; Bernoulli mean expected matches `PM·d`; default `PM = .80`, `PI = .10`. Implemented as `ComputeBernoulliStatistics`.

**Intentionally simplified:**

- The default `FindMicrosatellites` uses exact motif matching only; **consequence:** interrupted, impure, or mismatch-tolerant microsatellites are split into separate perfect tracts (use the opt-in `FindApproximateTandemRepeats` for those).
- Redundant-unit filtering and maximal-run reporting; **consequence:** each locus is reported once per primitive unit length (run-start phase), not once per rotation; compound / cross-size merging (MISA `interruptions`) is not performed.
- (Approximate, TRF [6][7]) **Partial detection criteria:** TRF's apparent-size / waiting-time criterion (cut-offs estimated by simulation), the random-walk distance range `d ± ⌊2.3·√(PI·d)⌋` summation, the best-period list for `d > 250` and the narrow-band WDP for patterns > 20 are not reproduced (a line-by-line port is also excluded by licence: TRF is AGPL-3.0, this library MIT). **Consequence (measured on 700 random sequences with embedded imperfect repeats, N and substitutions/indels, TRF `2 7 7 80 10 50 500`):** embedded periods ≤ 20 — 1 069 / 1 154 TRF rows (92.6 %) reproduced exactly (indices, period, consensus size, score, consensus), 96.0 % at region level (same period ±1, ≥ 90 % mutual overlap), C# rows confirmed by TRF 96.3 %; embedded periods ≤ 100 — 80.5 % exact, 93.1 % region level, 97.3 % confirmed. On exactly reproduced rows the numeric fields (copies, % matches, % indels, composition, entropy) agree on 1 347 / 1 349 (two equal-score traceback ties reached from a different trigger position). The gap comes from TRF's banded alignment and range distances (large periods) and trigger-position-dependent consensus choice. **Users should rely on:** compiled TRF [7] when bit-identical TRF output is required.
- (Approximate) `Entropy` uses the canonical ACGT-normalised Shannon entropy; TRF divides the four counts by the region length including N, so for regions containing N the TRF value is lower (TRF's is not a normalised entropy).

**Not implemented:**

- TRF's apparent-size criterion, random-walk distance ranges, best-period list and narrow-band WDP (see "Intentionally simplified"); TRF's alignment / flanking HTML outputs and masked-sequence file.
- PCR-stutter modeling and locus-specific forensic interpretation; **users should rely on:** dedicated forensic STR pipelines.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | No global non-overlap rule across unit lengths; same-period runs may overlap by < p. | Deviation | Overlapping loci of different unit lengths are all reported (MISA per-size convention); rotations of one run are never reported (fixed 2026-09). | accepted | The legacy doc described non-overlap broadly, but the source suppresses only contained intervals. |
| 2 | Approximate detector implements TRF's k-tuple + sum-of-heads trigger but not the apparent-size, range-distance, best-period-list criteria or the narrow band. | Deviation | Reported loci can differ from TRF (measured agreement in §5.3, REP-APPROX-001 Evidence); the analysis of a candidate is TRF-identical for patterns ≤ 20. | accepted (review 2026-09) | Replaces the former exhaustive O(n²·P·L²) window scan. |
| 3 | Percentages are exact (TRF truncates to integers); `Start` 0-based (TRF 1-based); output ordered by start (TRF: discovery order). | Convention | None numerically. | accepted | REP-APPROX-001 Evidence. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns empty enumerable. | There is no region long enough to contain the required repeat copies. |
| No repeats found | Returns empty enumerable. | No candidate motif reaches `minRepeats`. |
| Entire sequence is one repeat tract | Returns a result covering that tract. | The extension loop continues until the repeated pattern stops. |
| `minRepeats = 2` | Accepted. | The implementation rejects only values below `2`. |
| `minRepeats < 2` | Throws `ArgumentOutOfRangeException`. | Explicit parameter validation enforces this floor. |
| `maxUnitLength < minUnitLength` | Throws `ArgumentOutOfRangeException`. | Explicit parameter validation enforces ordered bounds. |
| `null` `DnaSequence` | Throws `ArgumentNullException`. | Explicit null guard on `DnaSequence` overloads. |
| Trailing partial copy (`ATATATA`) | `AT×3@0` only; the partial `A` is not counted and `TA×3@1` is not reported. | Maximal run, complete copies (pytrf `repeat = length / p`). |
| Run of `N` / non-ACGT unit (string overload) | Not reported. | MISA `[acgt]` motifs; pytrf skips `N`. |
| Cancellable overloads with invalid bounds | Throw `ArgumentOutOfRangeException` eagerly (all four overloads). | Shared validation (the DnaSequence cancellable overload previously skipped it; `minUnitLength = 0` never terminated). |
| Approximate: empty / too-short sequence | Returns empty enumerable. | No window of two copies exists. |
| Approximate: perfect tract | `PercentMatches = 100`, `PercentIndels = 0`, exact period/copy number (TRF README test_seqs tables reproduced). | A perfect alignment has only match columns. |
| Approximate: tract scoring below `minScore` | Not reported. | Benson (1999) report threshold [6]. |
| Approximate: short interrupted tract without a k-run (`CACACATACACA`) | Not reported (TRF reports nothing). | Sum-of-heads criterion 5 not met [6]. |
| Approximate: all-N / N inside a repeat | All-N: nothing; N inside: a mismatch. | TRF similarity matrix [7]. |
| Approximate: lowercase input | Same result as uppercase. | TRF uppercases input [7]. |
| Approximate: `minPeriod < 1`, `maxPeriod < minPeriod`, `maxPeriod > 2000`, `minScore < 1` | Throws `ArgumentOutOfRangeException` eagerly (both overloads, also for empty input). | Explicit parameter validation; TRF limits [7]. |

### 6.2 Limitations

The default detector detects only exact consecutive repeats and does not model interruptions, motif degeneracy, or sequencing noise; the opt-in `FindApproximateTandemRepeats` closes that gap for substitutions and indels with the TRF model (partial detection criteria, see §5.3). Both are limited to motif/period lengths within the configured range, which defaults to the biological microsatellite window of 1-6 bp. Default output is also canonicalized by redundant-unit filtering and contained-interval suppression, so it is not a complete enumeration of every equivalent motif interpretation.

## 7. Examples and Related Material

### 7.2 Related Use Cases

Trinucleotide repeat expansions highlighted in the legacy documentation include:

| Disease | Gene | Repeat | Normal | Pathogenic |
|---------|------|--------|--------|------------|
| Huntington's disease | `HTT` | `CAG` | 6-35 | 36-250 |
| Fragile X syndrome | `FMR1` | `CGG` | 6-53 | 230+ |
| Friedreich's ataxia | `FXN` | `GAA` | 7-34 | 100+ |
| Myotonic dystrophy 1 | `DMPK` | `CTG` | 5-34 | 50+ |

Additional common uses include forensic DNA profiling, where tetra- and pentanucleotide STR markers are preferred, and CODIS-style locus panels for identity testing [1][2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_Microsatellite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Microsatellite_Tests.cs)
- Approximate-detector tests: [RepeatFinder_ApproximateTandemRepeats_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_ApproximateTandemRepeats_Tests.cs) — TRF-locked rows; covers `INV-06`..`INV-09`
- Approximate-detector evidence: [REP-APPROX-001-Evidence.md](../../Evidence/REP-APPROX-001-Evidence.md); test spec [REP-APPROX-001.md](../../../tests/TestSpecs/REP-APPROX-001.md)
- Test spec: [REP-STR-001.md](../../../tests/TestSpecs/REP-STR-001.md)
- Related property tests: [RepeatFinderProperties.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Properties/RepeatFinderProperties.cs)
- Related metamorphic tests: [MetamorphicTests.cs](../../../tests/SuffixTree/SuffixTree.Tests/Algorithms/MetamorphicTests.cs)

## 8. References

1. Wikipedia. 2026. Microsatellite. Wikipedia. https://en.wikipedia.org/wiki/Microsatellite
2. Wikipedia. 2026. Trinucleotide repeat disorder. Wikipedia. https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder
3. Richard GF, Kerrest A, Dujon B. 2008. Comparative genomics and molecular dynamics of DNA repeats in eukaryotes. Microbiology and Molecular Biology Reviews. 72(4):686-727.
4. Tóth G, Gáspári Z, Jurka J. 2000. Microsatellites in different eukaryotic genomes: survey and analysis. Genome Research. 10(7):967-981.
5. Brinkmann B, Klintschar M, Neuhuber F, Hühne J, Rolf B. 1998. Mutation rate in human microsatellites. American Journal of Human Genetics.
6. Benson G. 1999. Tandem repeats finder: a program to analyze DNA sequences. Nucleic Acids Research. 27(2):573-580. https://doi.org/10.1093/nar/27.2.573
7. Benson G, Hernandez Y, Gelfand Y, Rodriguez A. Tandem Repeats Finder 4.10.0 — README ("TRF Definitions", "How does Tandem Repeats Finder work?") and source (`tr30dat.c`, `trfclean.h`; AGPL-3.0). https://github.com/Benson-Genomics-Lab/TRF (commit 355c1f9, 2020-06-29)
8. Kolpakov R, Kucherov G. 1999. Finding maximal repetitions in a word in linear time. Proc. 40th IEEE FOCS, 596-604. https://doi.org/10.1109/SFFCS.1999.814634
9. Thiel T, Michalek W, Varshney RK, Graner A. 2003. Exploiting EST databases for the development and characterization of gene-derived SSR-markers in barley. Theor Appl Genet 106:411-422 (MISA; source `misa.pl` v1.0, mirror https://raw.githubusercontent.com/cfljam/SSR_marker_design/master/misa.pl).
10. Du L, Zhang C, Liu Q, Zhang X, Yue B. 2018. Krait: an ultrafast tool for genome-wide survey of microsatellites and primer design. Bioinformatics 34(4):681-683 (pytrf 1.5.0, PyPI, `src/str.c`).
