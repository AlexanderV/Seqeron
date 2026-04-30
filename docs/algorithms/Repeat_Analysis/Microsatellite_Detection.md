# Microsatellite Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-STR-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Microsatellite detection identifies short tandem repeats (STRs, also called microsatellites or simple sequence repeats) whose motif length is between 1 and 6 nucleotides [1][3][4]. The repository implements exact consecutive-repeat detection in `RepeatFinder.FindMicrosatellites`, classifies each hit by repeat-unit length, and exposes overloads for `DnaSequence`, raw strings, and cancellation-aware execution. The implementation also removes redundant compound motifs such as `ATAT` when they are just repetitions of a smaller motif, and it suppresses results fully contained inside already reported repeat intervals. Microsatellite biology matters clinically, forensically, and evolutionarily because repeat expansions drive many genetic disorders and STR polymorphism underpins DNA profiling [1][2][3].

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

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every result satisfies `RepeatCount >= minRepeats`. | Results are emitted only after counting the required number of consecutive copies. |
| INV-02 | `minUnitLength <= RepeatUnit.Length <= maxUnitLength`. | The outer search loop enumerates only those unit lengths. |
| INV-03 | `TotalLength = RepeatUnit.Length × RepeatCount`. | `MicrosatelliteResult.TotalLength` is constructed from those two fields. |
| INV-04 | `RepeatType` matches the reported unit length. | `ClassifyRepeatType` maps unit lengths 1 through 6 to the corresponding repeat class. |
| INV-05 | Fully contained repeats are suppressed once a containing interval has already been reported. | The implementation checks previously reported `(Start, End)` intervals before yielding a new hit. |

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

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | 0-based start position of the reported microsatellite. |
| `RepeatUnit` | `string` | Repeated motif. |
| `RepeatCount` | `int` | Number of consecutive copies of the motif. |
| `TotalLength` | `int` | Total length of the repeat tract in bases. |
| `RepeatType` | `RepeatType` | Unit-length classification from mono- through hexanucleotide. |

### 3.3 Preconditions and Validation

All overloads validate `minUnitLength`, `maxUnitLength`, and `minRepeats`, rejecting `minUnitLength < 1`, `maxUnitLength < minUnitLength`, and `minRepeats < 2`. `DnaSequence` overloads throw `ArgumentNullException` for `null` sequence input. String overloads yield no results for `null` or empty input and normalize non-empty input to uppercase before scanning. Reported positions are 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize raw-string input to uppercase when needed.
2. For each unit length from `minUnitLength` to `maxUnitLength`, enumerate starting positions where at least `minRepeats` copies could fit.
3. Extract the candidate motif and skip it if it is redundant, meaning it is itself composed of a smaller repeated subunit.
4. Count consecutive occurrences of the motif.
5. If the repeat count is large enough and the resulting interval is not fully contained within an already reported interval, emit a `MicrosatelliteResult` with the appropriate `RepeatType`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Microsatellite detection | `O(n × U × R)` | `O(k)` | `U` is the searched unit-length range, `R` is the average repeat count encountered while extending motifs, and `k` is the number of retained intervals/results. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindMicrosatellites(DnaSequence, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable overload with progress reporting.
- `RepeatFinder.FindMicrosatellites(string, int, int, int)`: Raw-string overload with uppercase normalization.
- `RepeatFinder.FindMicrosatellites(string, int, int, int, CancellationToken, IProgress<double>?)`: Cancellable raw-string overload.

### 5.2 Current Behavior

The implementation filters redundant motifs with `IsRedundantUnit`, so larger patterns such as `ATAT` or `CAGCAG` are not reported when they can be expressed as repetitions of smaller units. It tracks previously reported `(Start, End)` intervals and suppresses a new result only when that new interval is completely contained within an existing reported interval; this is narrower than blanket overlap suppression. The cancellable implementation checks `cancellationToken` periodically and reports progress based on processed candidate positions. Raw-string overloads normalize input with `ToUpperInvariant()`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact detection of 1-6 bp short tandem repeats (STRs) [1][3][4].
- Classification of repeats into mono-, di-, tri-, tetra-, penta-, and hexanucleotide categories [1][3].
- Consecutive-copy counting with explicit reporting of repeat unit, count, start position, and total length.

**Intentionally simplified:**

- Exact motif matching only; **consequence:** interrupted, impure, or mismatch-tolerant microsatellites are not reported.
- Redundant-unit filtering and contained-interval suppression; **consequence:** output favors a canonical representative interval rather than an exhaustive list of every possible equivalent unit decomposition.

**Not implemented:**

- Approximate STR scoring, PCR-stutter modeling, and locus-specific forensic interpretation; **users should rely on:** `RepeatFinder.FindMicrosatellites` only for exact motif-repeat discovery.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Containment suppression is narrower than a global non-overlap rule. | Deviation | Some partially overlapping candidate repeats may still be reported if neither interval fully contains the other. | accepted | The legacy doc described non-overlap broadly, but the source suppresses only contained intervals. |

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

### 6.2 Limitations

The algorithm detects only exact consecutive repeats and does not model interruptions, motif degeneracy, or sequencing noise. It is limited to unit lengths within the configured range, which defaults to the biological microsatellite window of 1-6 bp. Output is also canonicalized by redundant-unit filtering and contained-interval suppression, so it is not a complete enumeration of every equivalent motif interpretation.

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

- Tests: [RepeatFinder_Microsatellite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_Microsatellite_Tests.cs)
- Test spec: [REP-STR-001.md](../../../tests/TestSpecs/REP-STR-001.md)
- Related property tests: [RepeatFinderProperties.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Properties/RepeatFinderProperties.cs)
- Related metamorphic tests: [MetamorphicTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MetamorphicTests.cs)

## 8. References

1. Wikipedia. 2026. Microsatellite. Wikipedia. https://en.wikipedia.org/wiki/Microsatellite
2. Wikipedia. 2026. Trinucleotide repeat disorder. Wikipedia. https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder
3. Richard GF, Kerrest A, Dujon B. 2008. Comparative genomics and molecular dynamics of DNA repeats in eukaryotes. Microbiology and Molecular Biology Reviews. 72(4):686-727.
4. Tóth G, Gáspári Z, Jurka J. 2000. Microsatellites in different eukaryotic genomes: survey and analysis. Genome Research. 10(7):967-981.
5. Brinkmann B, Klintschar M, Neuhuber F, Hühne J, Rolf B. 1998. Mutation rate in human microsatellites. American Journal of Human Genetics.
