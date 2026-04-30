# FASTQ Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | FileIO |
| Test Unit ID | PARSE-FASTQ-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

FASTQ parsing reads nucleotide sequences together with per-base quality scores.[1][2][3] In this repository, `FastqParser` parses FASTQ content from strings, files, and readers; decodes Phred quality encodings; filters and trims reads; computes summary statistics; and provides paired-end and writing helpers. The implementation supports both Phred+33 and Phred+64 encodings and can auto-detect between them using a character-range heuristic, but that auto mode is not reliable across the full Phred+33 printable range. It is simplified relative to the broader FASTQ ecosystem because it favors tolerant parsing and utility helpers over strict malformed-record rejection.[1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

FASTQ stores a sequence and its per-base quality scores in a text record whose canonical structure is four logical lines: a header beginning with `@`, a sequence, a separator beginning with `+`, and a quality string of the same length as the sequence.[1][2] Quality values are Phred scores encoded as printable ASCII characters.[1][2]

### 2.2 Core Model

The canonical record layout is:[1][2]

```text
@<identifier> <optional description>
<sequence>
+
<quality string>
```

Phred quality is defined as:[1][2]

$$
Q = -10 \cdot \log_{10}(p)
$$

with inverse:

$$
p = 10^{-Q/10}
$$

The encoding schemes preserved from the current document are:[1][2]

| Format | Offset | ASCII Range | Typical Q Range | Usage |
|--------|--------|-------------|-----------------|-------|
| Phred+33 | 33 | `!` to `~` | `0..93` | Sanger, Illumina 1.8+, PacBio, Nanopore |
| Phred+64 | 64 | `@` to `~` | `0..62` | Legacy Illumina 1.3-1.7 |

Example Phred values from the current document are:[1][2]

| Q Score | Error Probability | Accuracy |
|---------|-------------------|----------|
| 10 | 10% | 90% |
| 20 | 1% | 99% |
| 30 | 0.1% | 99.9% |
| 40 | 0.01% | 99.99% |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | In valid FASTQ, quality-string length equals sequence length | FASTQ encodes one quality score per base.[1][2] |
| INV-02 | Phred+33 scores lie in `0..93` and Phred+64 scores lie in `0..62` | These are the printable ASCII ranges of the two encodings.[1][2] |
| INV-03 | `Q = -10 log10(p)` and `p = 10^{-Q/10}` define the Phred/error-probability relationship | Standard Phred scoring model.[1][2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `filePath` | `string` | required | Path passed to `ParseFile(...)` | Missing or empty paths yield no records |
| `content` | `string` | required | FASTQ text passed to `Parse(...)` | Null or empty input yields no records |
| `reader` | `TextReader` | required | Reader passed to `Parse(...)` | Parsed line by line |
| `encoding` | `FastqParser.QualityEncoding` | `Auto` | Quality encoding used during decoding | `Auto` uses the repository's heuristic detector |
| `minAverageQuality` | `double` | required | Average-quality threshold for `FilterByQuality(...)` | Compared against decoded Phred values |
| `minQuality` | `int` | `20` | End-trimming threshold for `TrimByQuality(...)` | Applied to decoded scores |
| `adapter` | `string` | required | Adapter sequence for `TrimAdapter(...)` | No trimming occurs when adapter is null, empty, or shorter than `minOverlap` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Id` | `string` | Header identifier before the first space |
| `Description` | `string` | Remaining header text after the first space |
| `Sequence` | `string` | Parsed nucleotide sequence |
| `QualityString` | `string` | Parsed encoded quality string |
| `QualityScores` | `IReadOnlyList<int>` | Decoded Phred scores |
| `FastqStatistics.TotalReads` | `int` | Number of parsed records |
| `FastqStatistics.TotalBases` | `long` | Total bases across all records |
| `FastqStatistics.MeanReadLength` | `double` | Average sequence length |
| `FastqStatistics.MeanQuality` | `double` | Average decoded quality per base |
| `FastqStatistics.Q20Percentage` / `Q30Percentage` | `double` | Percent of bases with Q20/Q30 or better |
| `FastqStatistics.GcContent` | `double` | Fraction of bases that are `G` or `C` in the implementation |

### 3.3 Preconditions and Validation

The parser skips blank lines and only starts a record when a line begins with `@`. Sequence lines are accumulated until a line beginning with `+` is encountered. Quality lines are then accumulated until the quality string length reaches the sequence length. Null or empty input returns no records. `DecodeQualityScores(...)` returns an empty array for null or empty quality strings. `EncodeQualityScores(...)` clamps scores to the representable range of the selected encoding. `ErrorProbabilityToPhred(...)` returns `93` when the probability is zero or negative.

## 4. Algorithm

### 4.1 High-Level Steps

1. Read input until a header line beginning with `@` is found.
2. Split the header into `Id` and optional `Description` at the first space.
3. Accumulate sequence lines until the `+` separator line.
4. Accumulate quality lines until their total length reaches the sequence length.
5. Detect or apply the selected quality encoding.
6. Decode the quality string to Phred scores and yield a `FastqRecord`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository's auto-detection heuristic is:

```text
If any quality character is below '@', choose Phred+33.
Else if any quality character is above 'I', choose Phred+64.
Else default to Phred+33.
```

Paired-end support is modeled in two forms from the current document:[1]

| Mode | Description |
|------|-------------|
| Separate files | Read 1 and read 2 stored in separate files with matched order |
| Interleaved | Read 1 and read 2 alternate within a single FASTQ stream |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Parse` / `ParseFile` | `O(n)` | `O(1)` auxiliary | Linear scan over FASTQ text |
| `DetectEncoding` | `O(m)` | `O(1)` | `m` = quality-string length |
| `DecodeQualityScores` / `EncodeQualityScores` | `O(m)` | `O(m)` | Per-record quality conversion |
| `FilterByQuality` | `O(r * m)` | `O(1)` auxiliary | `r` = record count |
| `TrimByQuality` | `O(m)` | `O(m)` | End trimming over one record |
| `TrimAdapter` | `O(m * a)` worst case | `O(m)` | `a` = adapter length |
| `CalculateStatistics` | `O(r * m)` | `O(1)` auxiliary | Record-level aggregation |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [FastqParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastqParser.cs)

- `FastqParser.ParseFile(string, QualityEncoding)` and `Parse(string, QualityEncoding)`: Parse FASTQ records from files or text.
- `FastqParser.Parse(TextReader, QualityEncoding)`: Reader-based FASTQ parser.
- `FastqParser.DetectEncoding(...)`, `DecodeQualityScores(...)`, `EncodeQualityScores(...)`: Quality-encoding helpers.
- `FastqParser.PhredToErrorProbability(...)`, `ErrorProbabilityToPhred(...)`: Phred mathematics helpers.
- `FastqParser.FilterByQuality(...)`, `FilterByLength(...)`, `TrimByQuality(...)`, `TrimAdapter(...)`: Filtering and trimming utilities.
- `FastqParser.CalculateStatistics(...)`, `CalculatePositionQuality(...)`: Summary-statistics helpers.
- `FastqParser.WriteToStream(...)`, `ToFastqString(...)`, `InterleavePairedReads(...)`, `SplitInterleavedReads(...)`: Writing and paired-end helpers.

### 5.2 Current Behavior

`Parse(...)` is tolerant: it skips non-header lines, accepts multi-line sequences until `+`, and then reads quality lines until the accumulated quality length reaches the sequence length. It does not perform a strict malformed-record rejection pass beyond those rules. Auto-detection defaults ambiguous `@`-through-`I` ranges to Phred+33 and switches to Phred+64 whenever any quality character is above `I`, which means high-quality Phred+33 strings containing `J` through `~` are not distinguishable from Phred+64 by this heuristic alone. `DecodeQualityScores(...)` subtracts the selected ASCII offset and clamps negative values to `0`. `TrimByQuality(...)` trims only low-quality ends, not internal low-quality segments. `TrimAdapter(...)` first looks for an adapter overlap at the sequence end, then for a full adapter match that starts after position `0`; a full adapter already at the first base is left unchanged unless the end-overlap path trims it. `CalculateStatistics(...)` reports `Q20Percentage` and `Q30Percentage` as percentages, but `GcContent` as a `0..1` fraction.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- FASTQ parsing around `@` headers, `+` separator lines, and per-base quality strings.[1][2]
- Phred+33 and Phred+64 score conversion.[1][2]
- Phred/error-probability conversion formulas.[1][2]

**Intentionally simplified:**

- Auto encoding detection uses a character-range heuristic; **consequence:** ambiguous `@`-through-`I` ranges default to Phred+33, and high-quality Phred+33 strings containing `J` through `~` can be misclassified as Phred+64 unless callers specify the encoding explicitly.
- Parsing favors tolerant record assembly over strict format validation; **consequence:** malformed records can be skipped or partially assembled rather than rejected with an error.

**Not implemented:**

- Built-in support for FASTQ quality encodings beyond Phred+33 and Phred+64; **users should rely on:** no current alternative in this repository.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null input | Returns no records | Explicit early-return guards |
| Empty quality string | Detects as Phred+33 and decodes to an empty score list | Quality helpers guard empty input |
| All qualities below trim threshold | `TrimByQuality(...)` returns an empty-sequence record | End trimming can remove the full record |
| No adapter match | `TrimAdapter(...)` returns the original record | Adapter trimming is conditional |
| Interleaved input with odd record count | Final unmatched record stays on the alternating side reached by the splitter | `SplitInterleavedReads(...)` alternates records without pair validation |

### 6.2 Limitations

The repository provides practical FASTQ parsing and quality utilities, but it does not implement a strict FASTQ validator. Encoding detection is heuristic and cannot reliably disambiguate the full overlap between Phred+33 and Phred+64 printable ranges, and `GcContent` in summary statistics is a fraction rather than a percentage.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [FastqParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/FastqParserTests.cs)
- Related tests: [FastqParseTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqParseTests.cs), [FastqFilterTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqFilterTests.cs), [FastqStatisticsTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqStatisticsTests.cs), [FastqUtilityTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqUtilityTests.cs)
- Test specification: [PARSE-FASTQ-001.md](../../../tests/TestSpecs/PARSE-FASTQ-001.md)

## 8. References

1. Wikipedia contributors. FASTQ format. Wikipedia. https://en.wikipedia.org/wiki/FASTQ_format
2. Cock, P.J.A., et al. 2009. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research. https://doi.org/10.1093/nar/gkp1137
3. NCBI Sequence Read Archive. FASTQ and related submit formats. https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
