# Microsatellite Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-STR-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production (perfect STRs incl. MISA per-unit-size thresholds and compound SSRs, misa.pl-verified; opt-in approximate detector: TRF analysis and detection pipeline, 99.8–100 % of compiled-TRF rows identical — §5.3) |
| Last Reviewed | 2026-10-01 (B04 WP7) |

## 1. Overview

Microsatellite detection identifies short tandem repeats (STRs, also called microsatellites or simple sequence repeats) whose motif length is between 1 and 6 nucleotides [1][3][4]. The repository implements exact consecutive-repeat detection in `RepeatFinder.FindMicrosatellites` (the default), classifies each hit by repeat-unit length, and exposes overloads for `DnaSequence`, raw strings, and cancellation-aware execution. The implementation also removes redundant compound motifs such as `ATAT` when they are just repetitions of a smaller motif, and it reports each maximal primitive run once per unit length. MISA-style per-unit-size thresholds (`MisaDefaultMinRepeats` = `1-10 2-6 3-5 4-5 5-5 6-5`) and MISA compound microsatellites (SSRs interrupted by ≤ 100 bases, types `c` / `c*`) are available as additive overloads [9]. An **opt-in approximate detector**, `RepeatFinder.FindApproximateTandemRepeats`, additionally finds **imperfect / interrupted** tandem repeats (those containing substitutions or indels) using the Tandem Repeats Finder (TRF) model [6][7]: k-tuple matches at a common distance that pass Benson's sum-of-heads criterion trigger a wraparound-dynamic-programming (WDP) alignment of the sequence against tandem copies of the candidate pattern, the consensus is taken by majority rule, the sequence is realigned against it, and the TRF table (indices, period, copy number, consensus size, % matches and % indels between adjacent copies, score, composition, entropy) is reported (review 2026-09, REP-APPROX-001: identical to compiled TRF 4.10.0 on every analysed candidate with pattern ≤ 20). The default perfect-repeat detector is unchanged. Microsatellite biology matters clinically, forensically, and evolutionarily because repeat expansions drive many genetic disorders and STR polymorphism underpins DNA profiling [1][2][3].

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
| `progress` | `IProgress<double>?` | optional | Progress callback for cancellable scans. | Non-decreasing values in `[0, 1)` every 1000 visited run starts (fraction of the unit-length × position scan space), then exactly `1.0` when the scan completes; reported only while the result is enumerated (lazy). |
| `minRepeatsByUnitLength` | `IReadOnlyDictionary<int,int>` | — | Per-unit-size minimum copies (MISA `definition`), e.g. `RepeatFinder.MisaDefaultMinRepeats`; only the listed unit lengths are searched. | Non-null, non-empty (else `ArgumentNullException` / `ArgumentException`); unit lengths ≥ 1, copies ≥ 2 (else `ArgumentOutOfRangeException`); eager. |
| `maxInterruption` | `int` | `100` | (compound) Maximal number of bases between two consecutive SSRs of a compound (MISA `interruptions`, `MisaDefaultMaxInterruption`). | Negative values throw `ArgumentOutOfRangeException`. |
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
| `Entropy` | `double` | Shannon entropy (bits, 0–2) of the region's A/C/G/T composition, via canonical `SequenceComplexity.CalculateShannonEntropy` [7]; equals TRF's column unless the region contains N. |
| `EntropyTrf` | `double` | TRF's "Entropy (0-2)" column exactly: −Σ p_b log₂ p_b over A/C/G/T with p_b = count_b / `SpanLength` (N in the denominator, not in the sum) [7]. |
| `AlignedSequence` / `AlignedConsensus` | `string?` | The final alignment's sequence and consensus rows, left to right, '-' = gap (TRF alignment file) [7]. |
| `LeftFlank` / `RightFlank` | `string?` | Up to `FlankLength` symbols on each side (TRF `-f`: 500; `-ngs`: 50); null unless requested [7]. |

`TandemRepeatsFinderParameters` (TRF command line `Match Mismatch Delta PM PI Minscore MaxPeriod [-l n] [-r] [-f]`, defaults = recommended `2 7 7 80 10 50 500` [7]): `MatchWeight` (≥ 1), `MismatchPenalty` / `IndelPenalty` (≥ 1, applied negative; 3–7 recommended), `MatchProbability` (PM, 80 or 75 — the values TRF has data for), `IndelProbability` (PI 1..100; data for 10 / 20), `MinScore` (≥ 1), `MaxPeriod` (1..2000), `MaxRepeatLength` (`-l`, bp, default 2 000 000), `EliminateRedundancy` (`-r` = false), `FlankLength` (≥ 0). `MaskApproximateTandemRepeats(sequence, parameters?, softMask)` returns TRF's `-m` masked sequence (repeat positions → `N`, or lower case).

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

1. Scan positions `i` (1-based internally, as TRF) and distances `d = 1..D` (legacy overloads: `D = maxPeriod`; parameter-set API: TRF's MAXDISTANCE = max(200, min(max(MaxPeriod, 500), ⌊0.6·n⌋)), MaxPeriod then only filters the output). A *k-tuple match* at distance `d` ends at `i` when the last `k` symbols equal those `d` earlier (A/C/G/T only); PM = 80: `k = 4 / 5 / 7` for `d ≤ 29 / 30..159 / ≥ 160` (Benson 1999 Table 1); PM = 75: `k = 3 / 4 / 5 / 7` for `d ≤ 29 / 30..43 / 44..159 / ≥ 160` (TRF 4.10.0). Runs of adjacent tuple matches are kept in a window of the last `max(d, 20)` positions.
2. **Detection criteria** (Benson 1999; TRF README "Detection component"): (a) **sum of heads** — the heads in the window must reach `max(k+1, 5, ⌊μ − 1.65σ⌋)`, `μ, σ` = exact mean / s.d. of `R(d,k,PM)` (heads in runs ≥ k of a Bernoulli(PM) sequence of length `d`) — reproduces TRF's `sumdata80` and `sumdata75` tables for all `d ≤ 2000`; (b) **apparent size** — the first tuple match in the window must end at most `W = max(d,20) − y − 1` positions after the window's left end, where `y` is the largest number with `P(S > y | R ≥ x) ≥ 0.95` for `S` = distance between the first and last k-run of a Bernoulli(PM) sequence of length `max(d,20)` (TRF estimates this distribution by simulation; it is computed exactly here by a run-length Markov chain: README example PM .75 / k 5 / d 100 → 56). A distance `d > 20` that fails alone still qualifies through the **random-walk range** `d ± ⌊2.3·√(PI·d)⌋`: no *active* distance in the range may have more heads than `d`, and the heads summed over `d` and its lower range, or over a window of the same width sliding up through the upper range, must reach the cut-off while one summed distance with ≥ 35 % of `min(cut-off, heads(d))` passes the apparent-size test. A distance becomes active when it reaches a criteria test itself and inactive when a range scan finds its window empty. Positions already covered by an alignment at the same `d` are skipped. (c) For `d > 250`, the **best-period list**: if a region analysed earlier spans `i − 2d + 1 + W .. i`, `d` must be among that region's five best periods.
3. **WDP** with pattern `S[i−d+1..i]`: a backward local scan finds the leftmost row reaching the best score; a forward local alignment starting there gives the optimum. Patterns ≤ 20: full WDP (zero cells beyond `i` / before `i − max(d,20)` are killed; optimum = first strictly greatest cell). Patterns > 20: **narrow-band WDP** — each row keeps `2w+1` cells around a band centre that advances one pattern column per row and is moved onto the row maximum after 3 consecutive diagonal matches; `w = max(6, Δd_max)` backward, `min(2·max(6, Δd_max), ⌊size/3⌋)` forward (anchored on the backward optimum), kill rows `≤ i − size` / `≥ i`, optimum = last cell reaching the maximum, traceback continues through zero cells that are genuine path continuations. Traceback preference: match/mismatch, then sequence-vs-gap, then pattern-vs-gap.
4. Require ≥ 1.9 copies (ramp to 1.8 for 50–100, 1.8 above 100) and `d` among the **three best periods** of the aligned region (dinucleotide-distance histogram with its least-squares trend removed; period 1: ≥ 80 % one base); the region and its five best periods enter the best-period list.
5. **Consensus** by majority rule (a position is deleted when gaps are at least as frequent; an insertion point receives its most frequent base when insertions occur in ≥ half the passes), then **realign** against the consensus (step 3) and re-apply the copy rule and `minScore`.
6. **Statistics** from the final alignment: two cursors one consensus period apart compare each copy with the next (match / mismatch / indel), period = most common matching distance, composition and entropy over the region.
7. Drop periods above `maxPeriod`, sort by start, apply TRF **redundancy elimination** (unless `EliminateRedundancy = false`, TRF `-r`) (≥ 90 % overlap: same period with no higher score, or a multiple period with score ≤ 1.1×, is removed), drop periods below `minPeriod`, order by (start, end, period).

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
| Approximate detection (TRF) | `O(n × P + Σ_candidates (L + d) × min(d, b))` | `O(P² + L × min(c, b))` | `P` = examined distances; each candidate costs one WDP over its region `L` — full width for patterns ≤ 20 (extent-only pass, traceback matrix only for candidates passing the copy / best-period tests), band `b = 2w+1` for larger patterns (`c` = consensus size). Apparent-size + sum-of-heads tables: ≈ 0.15 s once per PM. 1 Mb, recommended set (WP7): 5.1 s (WP6: 15.1 s; compiled TRF 1.8 s). Measured (Release, this container): 100 kb random, `maxPeriod` 500 → 0.34 s; 100 kb of mixed embedded repeats → 1.3 s (545 rows; TRF 549 rows); pathological 100 kb perfect `(CA)n` with `maxPeriod` 500 → 54 s (compiled TRF 74 s — every even `d` aligns the whole array in both). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable overload with progress reporting.
- `RepeatFinder.FindMicrosatellites(string, int, int, int)`: Raw-string overload with uppercase normalization.
- `RepeatFinder.FindMicrosatellites(string, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable raw-string overload.
- `RepeatFinder.FindMicrosatellites(DnaSequence | string, IReadOnlyDictionary<int,int> minRepeatsByUnitLength, CancellationToken = default, IProgress<double>? = null)`: per-unit-size thresholds (MISA `definition`; `MisaDefaultMinRepeats`).
- `RepeatFinder.FindCompoundMicrosatellites(DnaSequence | string, IReadOnlyDictionary<int,int>? = null (MISA defaults), int maxInterruption = 100)` and `(…, int minRepeats, int maxInterruption = 100)`: MISA compound SSRs → `CompoundMicrosatelliteResult(Start, End, Components, Interruptions, IsOverlapping, Notation)` with `MisaType` `c`/`c*` and `Length`.
- `RepeatFinder.AssembleCompoundMicrosatellites(string, IEnumerable<MicrosatelliteResult>, int maxInterruption = 100)`: MISA's compound chaining applied to any SSR list (stable by start; ties keep input order).
- `RepeatFinder.FindApproximateTandemRepeats(DnaSequence | string, int minPeriod, int maxPeriod, int minScore)`: Opt-in approximate / imperfect / interrupted tandem-repeat detector (TRF model [6][7]); wraparound DP is its own algorithm (the library's pairwise aligners have no wraparound mode); entropy via canonical `SequenceComplexity.CalculateShannonEntropy`.
- `RepeatFinder.FindApproximateTandemRepeats(DnaSequence | string, TandemRepeatsFinderParameters, int minPeriod = 1)`: the same detector with the full TRF parameter set (weights, PM 80/75, PI random-walk range, Minscore, MaxPeriod with TRF's MAXDISTANCE, `-l`, `-r`, `-f` flanks).
- `RepeatFinder.MaskApproximateTandemRepeats(string, TandemRepeatsFinderParameters? = null, bool softMask = false)` / `(string, IEnumerable<ApproximateTandemRepeatResult>, bool softMask = false)`: TRF `-m` masked sequence (repeat positions → `N`) or soft mask (lower case).
- `RepeatFinder.ComputeBernoulliStatistics(string repeatTract, int period, double expectedMatchProbability = 0.80)`: TRF Bernoulli measures (Benson 1999 [6]) — PM, PI and `ExpectedMatches = PM × trials` (= matches), computed **between adjacent copies** through the same TRF WDP + consensus alignment of the tract; flags whether PM ≥ the assessed PM (default 0.80).
- `RepeatFinder.TrfSumOfHeadsCriterion(int d[, int pm])` (internal): Benson's sum-of-heads cut-off from the exact moments of `R(d,k,PM)` (PM 80 or 75).

### 5.2 Current Behavior

Per unit length, each maximal perfect run is reported once at its left end with the run-start motif phase and the number of complete copies (review 2026-09, REP-STR-001). Redundant (non-primitive) units such as `ATAT` or `CAGCAG` are not reported — their runs appear at the primitive unit length. Units containing non-ACGT symbols (e.g. `N` in the raw-string overload) are never reported. Runs of different unit lengths may overlap and are reported independently (MISA per-motif-size convention); two distinct runs of the same period may overlap by fewer than `p` bases (e.g. `ACACACGCGCGC` → `AC×3@0`, `CG×3@5`), where greedy tools (MISA regex `/g`, pytrf `next_start`) instead restart after the first run's full copies (`GC×3@6`).

**Compound SSRs (MISA `misa.pl` v1.0 [9]).** SSRs are ordered by start; SSR i+1 joins the compound of SSR i when `start(i+1) − end(i) ≤ maxInterruption` (0-based start, exclusive end; adjacent = 0 bases and overlapping < 0 always join). The comparison uses the previous SSR's end, not the running maximum, and the compound's end is the last component's end (MISA `end` column; smaller than the span when the last component is nested). Type `c*` when any consecutive pair overlaps, else `c`; notation `(M)r` per component, the lower-cased interrupting bases between them, `*` after an overlapping component (e.g. `(TA)6tccgt(GA)7ttttt(A)12`, `(AC)7(CAG)6*`). **Cross-check (2026-09-30):** 6 048 sequences in 6 `misa.ini` configurations (`1-10 2-6 3-5 4-5 5-5 6-5` with interruptions 100 and 0; `1-5 2-3 …` 20; `1-3 2-2 …` 5; `2-4 3-3 5-2` 50; defaults + `7-4 8-3 10-3` 100), 72 974 misa.pl SSRs, 12 330 compounds: (a) this library's SSR lists = brute-force maximal primitive runs with per-size thresholds, 0 mismatching sequences; (b) `AssembleCompoundMicrosatellites` fed with misa.pl's own SSR list in misa.pl order reproduces every misa.pl `.misa` row (type, notation, size, start, end), 0 mismatching sequences; (c) end to end (this library's SSRs), `.misa` rows are identical in every sequence whose SSR list equals misa.pl's; the SSR-list differences (97 / 60 / 179 / 580 / 31 / 24 sequences per configuration) are all one of the two documented detection conventions — same-period run overlapping the previous match by < p bases (MISA truncates it: 1 220 SSRs) or a match of a non-primitive unit consuming bases before being rejected (1 365 SSRs) — 0 unexplained; plus 53 sequences (lowest-threshold configurations only) whose SSR sets agree but contain equal-start SSRs, which misa.pl orders by Perl hash order.

**Cross-check (2026-09):** 3,132 random + crafted cases (unit lengths 1–6, `minRepeats` 2–5, alphabets incl. `N`) agree 3,132/3,132 with an independent brute-force maximal-repetition reference; MISA-regex agreement 3,099/3,132 and pytrf 1.5.0 `STRFinder` (single motif size) 2,752/3,132, every disagreement being one of the two documented conventions (same-period runs overlapping by < p: 17; greedy tools consuming a non-primitive/`N` region first: MISA 16, pytrf 363 — pytrf has no primitivity check).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of 1-6 bp short tandem repeats (STRs) [1][3][4].
- Classification of repeats into mono-, di-, tri-, tetra-, penta-, and hexanucleotide categories [1][3].
- Consecutive-copy counting with explicit reporting of repeat unit, count, start position, and total length.

- (Approximate, TRF [6][7]) Analysis component: WDP local alignment, majority consensus, realignment, TRF statistics (indices, period, copies, consensus size, adjacent-copy % matches / % indels, score, composition), minimum copy number, three-best-periods test, redundancy elimination. **Cross-check (2026-09):** on 2 483 TRF-detected candidates (instrumented TRF 4.10.0, same position and distance) the C# analysis is identical (indices, score, consensus size, copy number, adjacent-copy counts) for all 1 524 candidates with pattern ≤ 20 (TRF's full-WDP range) and for 827 / 959 larger ones (TRF aligns those in a narrow band — implemented in WP7: alignment rows for consensus > 20 now 149/149 identical to TRF's alignment file).
- (Approximate, TRF [6]) Detection: k-tuple trigger with Benson's tuple sizes and the sum-of-heads criterion (exact moments + normal approximation, = TRF table for all d ≤ 2000).
- (Approximate, TRF [6]) Recommended scoring constants match `+2`, mismatch `−7`, indel `−7`, and the `Minscore = 50` report threshold; every TRF parameter is settable through `TandemRepeatsFinderParameters` (weights, PM 80/75, PI, Minscore, MaxPeriod, `-l`, `-r`, `-f`), the masked sequence (`-m`) through `MaskApproximateTandemRepeats`. **Cross-check (2026-10, WP7, 700 sequences):** 99.8–100 % of TRF rows exact / 100 % region level on all seven parameter sets (WP6: 87.8 % / 99.4 % recommended, 64.2 % / 93.8 % for 2 3 3 80 20); `-r` 99.8–99.9 %, `-l` 60/120/250 bp 100/100/99.8 %; masks identical for 699–700/700 sequences; `EntropyTrf`, 50/500-bp flanks and alignment rows (192/192 consensus ≤ 20, 149/149 above) identical to TRF on every same-locus row (REP-APPROX-001 Evidence, WP6/WP7 revisions).

- (Bernoulli, TRF [6]) The probabilistic measures: "We model alignment of two tandem copies … by … independent Bernoulli trials"; PM = P(Heads) = "the average percent identity between the copies"; PI = "the average percentage of insertions and deletions between the copies"; statistics "between adjacent copies … not between the sequence and the consensus pattern"; Bernoulli mean expected matches `PM·d`; default `PM = .80`, `PI = .10`. Implemented as `ComputeBernoulliStatistics`.

**Intentionally simplified:**

- The default `FindMicrosatellites` uses exact motif matching only; **consequence:** interrupted, impure, or mismatch-tolerant microsatellites are split into separate perfect tracts (use the opt-in `FindApproximateTandemRepeats` for those).
- Redundant-unit filtering and maximal-run reporting; **consequence:** each locus is reported once per primitive unit length (run-start phase), not once per rotation. Compound SSRs (MISA `interruptions`) are a separate, opt-in result (`FindCompoundMicrosatellites`), not a merge of the SSR list.
- (Approximate, TRF [6][7]) **Apparent-size cut-offs are exact, TRF's are simulated:** TRF ships Monte-Carlo estimates of the apparent-size criterion (`waitdata80/75`); this library computes the same distribution exactly (equal at 825/2000 and 713/2000 distances, |Δ| ≤ 3 / 4 elsewhere). **Consequence (WP7, 700 random sequences with embedded imperfect repeats, seven TRF parameter sets):** 99.8–100 % of TRF rows identical in every field, 100 % at region level; each of the remaining 16 rows traces to a single table entry where TRF's noisy value is 1 below the exact one (with TRF's table substituted: 100 %). A line-by-line port is excluded by licence (TRF AGPL-3.0, this library MIT); the implementation follows the published method and the README. TRF exits when a band exceeds 150 cells (PI 20, patterns ≥ 1365); no such limit here.
- (Approximate) `Entropy` uses the canonical ACGT-normalised Shannon entropy; TRF divides the four counts by the region length including N, so for regions containing N the TRF value is lower (TRF's is not a normalised entropy). `EntropyTrf` reports TRF's value in every case (WP6).

**Not implemented:**

- TRF's HTML formatting (the alignment rows, flanks and masked sequence themselves are exposed).
- PCR-stutter modeling and locus-specific forensic interpretation; **users should rely on:** dedicated forensic STR pipelines.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | No global non-overlap rule across unit lengths; same-period runs may overlap by < p. | Deviation | Overlapping loci of different unit lengths are all reported (MISA per-size convention); rotations of one run are never reported (fixed 2026-09). | accepted | The legacy doc described non-overlap broadly, but the source suppresses only contained intervals. |
| 4 | Detection convention vs misa.pl: maximal primitive runs instead of MISA's resumed regex scan. | Convention | misa.pl truncates a same-size run that overlaps the previous match by < p bases and hides a primitive run that starts inside a rejected non-primitive match; this library reports the maximal runs (e.g. `(AAGAA)8` at 45 vs misa.pl `(AGAAA)8` at 46). Compound assembly from a given SSR list is misa.pl-identical. | accepted (documented, tested) | §5.2 cross-check. |
| 5 | Compound SSRs with equal start: stable order (input order; `FindCompoundMicrosatellites` → unit length). | Convention | misa.pl uses Perl hash order (non-deterministic across runs). | accepted | Feeding misa.pl's order reproduces its rows. |
| 2 | Approximate detector implements TRF's full detection pipeline (k-tuples, sum of heads, apparent size, random-walk range over active distances, best-period list) and the narrow band; apparent-size cut-offs are exact instead of TRF's simulated table. | Convention | 99.8–100 % of TRF rows identical (§5.3; REP-APPROX-001 Evidence §WP7); every remaining row traced to one simulated table entry. | resolved (WP7) | Replaces the former exhaustive O(n²·P·L²) window scan (2026-09) and the partial detector (WP6). |
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

The default detector detects only exact consecutive repeats and does not model interruptions, motif degeneracy, or sequencing noise; the opt-in `FindApproximateTandemRepeats` closes that gap for substitutions and indels with the TRF model (full TRF detection pipeline; measured parity in §5.3). Both are limited to motif/period lengths within the configured range, which defaults to the biological microsatellite window of 1-6 bp. Default output is also canonicalized by redundant-unit filtering and maximal-run reporting, so it is not a complete enumeration of every equivalent motif interpretation.

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

- Tests: [RepeatFinder_Microsatellite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Microsatellite_Tests.cs), [RepeatFinder_MisaCompound_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_MisaCompound_Tests.cs) (MISA thresholds / compound SSRs / classes / progress)
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
