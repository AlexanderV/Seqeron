# Microsatellite Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-STR-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-24 |

## 1. Overview

Microsatellite detection identifies short tandem repeats (STRs, also called microsatellites or simple sequence repeats) whose motif length is between 1 and 6 nucleotides [1][3][4]. The repository implements exact consecutive-repeat detection in `RepeatFinder.FindMicrosatellites` (the default), classifies each hit by repeat-unit length, and exposes overloads for `DnaSequence`, raw strings, and cancellation-aware execution. The implementation also removes redundant compound motifs such as `ATAT` when they are just repetitions of a smaller motif, and it suppresses results fully contained inside already reported repeat intervals. An **opt-in approximate detector**, `RepeatFinder.FindApproximateTandemRepeats`, additionally finds **imperfect / interrupted** tandem repeats (those containing substitutions or indels) using the Tandem Repeats Finder alignment model [6]: a candidate pattern of each period is aligned against tandem copies of itself, the consensus is taken by majority rule, and the resulting alignment yields the period size, copy number, percent matches, percent indels, consensus pattern, and alignment score. The default perfect-repeat detector is unchanged. Microsatellite biology matters clinically, forensically, and evolutionarily because repeat expansions drive many genetic disorders and STR polymorphism underpins DNA profiling [1][2][3].

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

**Approximate (TRF) model.** Benson (1999) defines a tandem repeat as "two or more contiguous, *approximate* copies of a pattern of nucleotides" and reports, for each repeat, the period size, the number of copies aligned with the consensus pattern, the consensus size, the percent of matches and percent of indels between adjacent copies overall, and an alignment score [6]. The alignment is scored Smith-Waterman style with weights for match, mismatch and indels; the recommended parameter set is match `+2`, mismatch `7`, indel (delta) `7` (applied as negatives), and only repeats scoring at least `Minscore = 50` are reported [6]. The consensus pattern is determined "by majority rule from the alignment" [6]. The opt-in detector reproduces these statistics by aligning the observed window against a whole number of tandem copies of the majority-rule consensus and reading the match / mismatch / indel columns of the resulting alignment.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every result satisfies `RepeatCount >= minRepeats`. | Results are emitted only after counting the required number of consecutive copies. |
| INV-02 | `minUnitLength <= RepeatUnit.Length <= maxUnitLength`. | The outer search loop enumerates only those unit lengths. |
| INV-03 | `TotalLength = RepeatUnit.Length × RepeatCount`. | `MicrosatelliteResult.TotalLength` is constructed from those two fields. |
| INV-04 | `RepeatType` matches the reported unit length. | `ClassifyRepeatType` maps unit lengths 1 through 6 to the corresponding repeat class. |
| INV-05 | Fully contained repeats are suppressed once a containing interval has already been reported. | The implementation checks previously reported `(Start, End)` intervals before yielding a new hit. |
| INV-06 | Every approximate result has `AlignmentScore >= minScore`. | A candidate window is retained only when its TRF alignment score reaches the threshold [6]. |
| INV-07 | For an approximate result, `0 <= PercentMatches <= 100` and `0 <= PercentIndels <= 100`, and a perfect tract yields `PercentMatches = 100`, `PercentIndels = 0`. | Each percentage is a column count divided by the total alignment-column count; a perfect alignment has only match columns. |
| INV-08 | For an approximate result, `CopyNumber = (non-gap aligned bases) / Period` and `ConsensusSize = Period`. | Copy number is the aligned observed length over the period [6]; the majority-rule consensus is built with exactly `Period` columns. |

| INV-09 | (Bernoulli) For `ComputeBernoulliStatistics`, `MatchProbability ∈ [0,1]`, `IndelProbability ∈ [0,1]`, and `Matches + Mismatches + Indels = BernoulliTrials`; a perfect tract yields `MatchProbability = 1`, `IndelProbability = 0`. | Each Bernoulli trial is one alignment column between two adjacent copies, classified as exactly match / mismatch / indel [6]. |
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
| `minPeriod` | `int` | `1` | (approximate) Minimum period (motif) size to consider. | `FindApproximateTandemRepeats`; values below `1` throw `ArgumentOutOfRangeException`. |
| `maxPeriod` | `int` | `6` | (approximate) Maximum period (motif) size to consider. | `FindApproximateTandemRepeats`; values below `minPeriod` throw `ArgumentOutOfRangeException`. |
| `minScore` | `int` | `50` | (approximate) Minimum TRF alignment score to report. | `FindApproximateTandemRepeats`; default `DefaultApproximateMinScore = 50` per Benson (1999) [6]. |
| `repeatTract` | `string` | required | (Bernoulli) Observed tandem-repeat tract (≥ 2 copies). | `ComputeBernoulliStatistics`; `null` throws `ArgumentNullException`; fewer than two copies throws `ArgumentException`. |
| `period` | `int` | required | (Bernoulli) Repeat period (copy length). | `ComputeBernoulliStatistics`; values below `1` throw `ArgumentOutOfRangeException`. |
| `expectedMatchProbability` | `double` | `0.80` | (Bernoulli) PM threshold the tract is assessed against. | `ComputeBernoulliStatistics`; values outside `[0,1]` throw `ArgumentOutOfRangeException`; default `TrfDefaultMatchProbability = 0.80` per Benson (1999) [6]. |

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
| `Start` | `int` | 0-based start position of the approximate repeat window. |
| `SpanLength` | `int` | Number of observed (non-gap) bases spanned by the repeat. |
| `Period` | `int` | Period (motif) size. |
| `ConsensusSize` | `int` | Size of the consensus pattern (equals `Period` in this subset). |
| `Consensus` | `string` | Majority-rule consensus motif [6]. |
| `CopyNumber` | `double` | Copies aligned with the consensus = aligned bases / period [6]. |
| `PercentMatches` | `double` | Percent of matches between adjacent copies overall (0–100) [6]. |
| `PercentIndels` | `double` | Percent of indels between adjacent copies overall (0–100) [6]. |
| `AlignmentScore` | `int` | TRF alignment score (sum of column weights) [6]. |

`ComputeBernoulliStatistics` returns `TandemRepeatBernoulliStatistics` (the TRF Bernoulli statistical measures, Benson 1999 [6]):

| Field | Type | Description |
|-------|------|-------------|
| `Period` | `int` | Repeat period (copy length). |
| `AdjacentCopyPairs` | `int` | Number of adjacent copy pairs compared. |
| `BernoulliTrials` | `int` | Total Bernoulli trials = total alignment columns over all adjacent pairs [6]. |
| `Matches` / `Mismatches` / `Indels` | `int` | Column counts (heads = matches) over all adjacent-copy alignments. |
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
2. For each unit length from `minUnitLength` to `maxUnitLength`, enumerate starting positions where at least `minRepeats` copies could fit.
3. Extract the candidate motif and skip it if it is redundant, meaning it is itself composed of a smaller repeated subunit.
4. Count consecutive occurrences of the motif.
5. If the repeat count is large enough and the resulting interval is not fully contained within an already reported interval, emit a `MicrosatelliteResult` with the appropriate `RepeatType`.

For the approximate detector (`FindApproximateTandemRepeats`):

1. For each period from `minPeriod` to `maxPeriod` and each start position, grow a tandem window one base at a time (at least two copies).
2. Build the majority-rule consensus over the period-aligned columns of the window [6].
3. Align the window against a whole number of tandem copies of the consensus using TRF scoring (match `+2`, mismatch `−7`, indel `−7`) [6].
4. Read the alignment columns to compute matches, mismatches, indels → percent matches, percent indels, copy number, and alignment score.
5. Keep the best-scoring window per (start, period) when its score reaches `minScore`; report best-score-first, suppressing windows contained in a higher-scoring accepted window.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Approximate-detector scoring constants (Tandem Repeats Finder recommended set "2 7 7 … 50 …" [6]):

| Constant | Value | Source |
|----------|-------|--------|
| Match weight | `+2` | Benson (1999): match weight "+2 in all options" [6] |
| Mismatch penalty | `−7` | Benson (1999): recommended Mismatch = 7 [6] |
| Indel (delta) penalty | `−7` per gap column | Benson (1999): recommended Delta = 7; flat per-column indel [6] |
| `DefaultApproximateMinScore` | `50` | Benson (1999): "Only those repeats scoring at least 50 … are reported" [6] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Microsatellite detection (perfect) | `O(n × U × R)` | `O(k)` | `U` is the searched unit-length range, `R` is the average repeat count encountered while extending motifs, and `k` is the number of retained intervals/results. |
| Approximate detection (TRF) | `O(n² × P × L²)` worst case | `O(L²)` | `P` is the period range and `L` the window length; each (start, period, window) does a Needleman-Wunsch alignment. Deterministic exhaustive scan; appropriate for short tracts, not whole-genome scale. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable overload with progress reporting.
- `RepeatFinder.FindMicrosatellites(string, int, int, int)`: Raw-string overload with uppercase normalization.
- `RepeatFinder.FindMicrosatellites(string, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable raw-string overload.
- `RepeatFinder.FindApproximateTandemRepeats(DnaSequence | string, int minPeriod, int maxPeriod, int minScore)`: Opt-in approximate / imperfect / interrupted tandem-repeat detector (TRF model [6]); reuses [`SequenceAligner.GlobalAlign`](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs) for the underlying alignment.
- `RepeatFinder.ComputeBernoulliStatistics(string repeatTract, int period, double expectedMatchProbability = 0.80)`: Opt-in TRF Bernoulli statistical-significance measures (Benson 1999 [6]) — PM (matching probability), PI (indel probability), and the Bernoulli-mean expected matches `PM·d`, estimated **between adjacent copies** (not against the consensus); flags whether the tract is at least as conserved as a random tandem repeat with the default PM = 0.80.

### 5.2 Current Behavior

The implementation filters redundant motifs with `IsRedundantUnit`, so larger patterns such as `ATAT` or `CAGCAG` are not reported when they can be expressed as repetitions of smaller units. It tracks previously reported `(Start, End)` intervals and suppresses a new result only when that new interval is completely contained within an existing reported interval; this is narrower than blanket overlap suppression. The cancellable implementation checks `cancellationToken` periodically and reports progress based on processed candidate positions. Raw-string overloads normalize input with `ToUpperInvariant()`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of 1-6 bp short tandem repeats (STRs) [1][3][4].
- Classification of repeats into mono-, di-, tri-, tetra-, penta-, and hexanucleotide categories [1][3].
- Consecutive-copy counting with explicit reporting of repeat unit, count, start position, and total length.

- (Approximate, TRF [6]) Reported statistics — period size, copy number, consensus size, percent matches and percent indels between adjacent copies overall, alignment score — computed from an alignment of the sequence against tandem copies of the majority-rule consensus.
- (Approximate, TRF [6]) Recommended scoring constants match `+2`, mismatch `−7`, indel `−7`, and the `Minscore = 50` report threshold.

- (Bernoulli, TRF [6]) The probabilistic measures: "We model alignment of two tandem copies … by … independent Bernoulli trials"; PM = P(Heads) = "the average percent identity between the copies"; PI = "the average percentage of insertions and deletions between the copies"; statistics "between adjacent copies … not between the sequence and the consensus pattern"; Bernoulli mean expected matches `PM·d`; default `PM = .80`, `PI = .10`. Implemented as `ComputeBernoulliStatistics`.

**Intentionally simplified:**

- The default `FindMicrosatellites` uses exact motif matching only; **consequence:** interrupted, impure, or mismatch-tolerant microsatellites are split into separate perfect tracts (use the opt-in `FindApproximateTandemRepeats` for those).
- Redundant-unit filtering and contained-interval suppression; **consequence:** output favors a canonical representative interval rather than an exhaustive list of every possible equivalent unit decomposition.
- (Approximate, TRF [6]) Candidate repeats are found by a **deterministic exhaustive (start, period) scan with alignment scoring**, in place of TRF's probabilistic k-tuple distance-list seeding; **consequence:** the reported statistics of a repeat are faithful to Benson (1999), but the candidate-discovery heuristic differs (the subset examines all windows up to `maxPeriod`, which limits practical sequence/period size rather than scaling to whole genomes).

- (Bernoulli, TRF [6]) PM/PI are estimated **between adjacent copies** by segmenting the tract into period-length copies and aligning each adjacent pair; **consequence:** for substitution/perfect tracts the estimate is exact, while for indel-containing tracts the per-pair alignment frame is alignment-dependent (the qualitative Bernoulli partition still holds, but the exact PI is frame-sensitive).

**Not implemented:**

- TRF's probabilistic k-tuple **seeding** — the sum-of-heads percentile cut-off `R(d,k,pM)` ("the largest x such that 95% of the time R(d,k,pM) ≥ x") and the random-walk band `W(d,pI)` — which drive whole-genome-scale candidate discovery. The per-repeat Bernoulli statistical measures (PM, PI, expected matches) ARE now computed by `ComputeBernoulliStatistics`; the residual is the **genome-scale-performance** seeding index (the deterministic exhaustive scan is not a seeded genome-scale index), whose 95% percentile cut-offs come from TRF's non-redistributable simulation tables; **users should rely on:** the reference Tandem Repeats Finder tool [6][7] for genome-scale seeded detection.
- PCR-stutter modeling and locus-specific forensic interpretation; **users should rely on:** dedicated forensic STR pipelines.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Containment suppression is narrower than a global non-overlap rule. | Deviation | Some partially overlapping candidate repeats may still be reported if neither interval fully contains the other. | accepted | The legacy doc described non-overlap broadly, but the source suppresses only contained intervals. |
| 2 | Approximate detector uses exhaustive (start, period) scanning, not TRF k-tuple seeding. | Assumption | Candidate discovery is deterministic but O(n²·P·L²); not whole-genome scale. Reported statistics are unaffected. | accepted | Honest residual; see §5.3 "Not implemented" and Evidence ASSUMPTION 1. Use the TRF tool [7] at genome scale. |
| 3 | Percent matches / percent indels use total alignment columns as the denominator. | Assumption | Reproduces Benson (1999) worked statistics; the source names the statistics but gives no verbatim percentage formula. | accepted | See Evidence ASSUMPTION 2. |

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
| Approximate: empty / too-short sequence | Returns empty enumerable. | No window of two copies exists. |
| Approximate: perfect tract | `PercentMatches = 100`, `PercentIndels = 0`, exact period/copy number. | A perfect alignment has only match columns. |
| Approximate: tract scoring below `minScore` | Not reported. | Benson (1999) report threshold [6]. |
| Approximate: `minPeriod < 1` or `maxPeriod < minPeriod` | Throws `ArgumentOutOfRangeException`. | Explicit parameter validation. |

### 6.2 Limitations

The default detector detects only exact consecutive repeats and does not model interruptions, motif degeneracy, or sequencing noise; the opt-in `FindApproximateTandemRepeats` closes that gap for substitutions and indels but uses an exhaustive period scan that is not whole-genome scale and does not compute TRF's statistical significance. Both are limited to motif/period lengths within the configured range, which defaults to the biological microsatellite window of 1-6 bp. Default output is also canonicalized by redundant-unit filtering and contained-interval suppression, so it is not a complete enumeration of every equivalent motif interpretation.

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
- Approximate-detector tests: [RepeatFinder_ApproximateTandemRepeats_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_ApproximateTandemRepeats_Tests.cs) — covers `INV-06`, `INV-07`, `INV-08`
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
7. Benson G. Tandem Repeats Finder — reference implementation and documentation. https://github.com/Benson-Genomics-Lab/TRF and https://tandem.bu.edu/trf/trf.definitions.html
