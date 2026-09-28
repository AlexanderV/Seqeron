# FASTQ Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | FileIO |
| Test Unit ID | PARSE-FASTQ-001 |
| Related Projects | N/A |
| Implementation Status | Complete (documented limitations) |
| Last Reviewed | 2026-09-28 |

## 1. Overview

FASTQ parsing reads nucleotide sequences together with per-base quality scores.[1][2][3] In this repository, `FastqParser` parses FASTQ content from strings, files, and readers; decodes Phred quality encodings; filters and trims reads; computes summary statistics; and provides paired-end and writing helpers. The implementation supports both Phred+33 and Phred+64 encodings; in auto mode the offset is detected once per file (FastQC-style, from the whole file's quality characters) by the canonical `QualityScoreAnalyzer.DetectEncoding(IEnumerable<string>)`. Record parsing follows Biopython's `FastqGeneralIterator` and rejects malformed records with `FormatException`.[1][2]

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

Blank lines between records are skipped; any other line where a record is expected must begin with `@`. Sequence lines are accumulated (right-trimmed) until the `+` line; whitespace inside the sequence is rejected. If the `+` line repeats a title it must equal the `@` title. Quality lines are accumulated by length — a line beginning with `@` only starts the next record once the quality is complete — and the quality length must equal the sequence length. Every violation, and truncation at end of file, throws `FormatException` (the same cases Biopython 1.88 `FastqGeneralIterator` rejects with `ValueError`). Null or empty input returns no records. `DecodeQualityScores(...)` returns an empty array for null or empty quality strings and throws `ArgumentOutOfRangeException` for a symbol outside the encoding's range (Biopython `InvalidCharError`); inside `Parse` this surfaces as `FormatException`. `EncodeQualityScores(...)` clamps scores to the representable range of the selected encoding (Biopython truncates at 93 / 62 likewise). `ErrorProbabilityToPhred(...)` returns `93` for p = 0, caps every result at `93`, and throws for p outside [0, 1] or NaN.

## 4. Algorithm

### 4.1 High-Level Steps

1. Skip blank lines; the next line must begin with `@` (else `FormatException`).
2. Split the title into `Id` and `Description` at the first whitespace character (Biopython `title.split(None, 1)`; shared `SequenceFormatHelper.SplitTitle`, also used by `FastaParser`).
3. Accumulate sequence lines until the `+` separator line; validate the optional `+` caption and the absence of whitespace.
4. Accumulate quality lines by length (an `@` line ends the record only once the quality is complete); require quality length = sequence length.
5. In Auto mode, determine the encoding once for the whole input (`ParseFile` / `Parse(string)`: a lightweight first pass; `Parse(TextReader)`: records are buffered because a reader cannot be rewound).
6. Decode the quality string to Phred scores (canonical `QualityScoreAnalyzer.ParseQualityString`) and yield a `FastqRecord`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The auto-detection rule (canonical `QualityScoreAnalyzer.DetectEncoding`, applied to the whole file in `Parse`/`ParseFile` and to one string in `DetectEncoding(string)`):

```text
If the lowest quality character is below ASCII 64 ('@'), choose Phred+33 (proven: outside Phred+64).
Else if the highest character is above ASCII 74 ('J' = Q41, Illumina 1.8+ Phred+33 ceiling), choose Phred+64 (inferred).
Else (all characters in ASCII 64-74) default to Phred+33 (ambiguous; LimitationPolicy "PARSE-FASTQ-001").
```

`TrimAdapter` implements cutadapt's regular 3' adapter removal with exact matching (`-a ADAPTER -e 0 -O minOverlap`): the read is cut at the leftmost position from which the read matches the adapter (a full occurrence anywhere, including position 0, or an adapter prefix of at least `minOverlap` bases running off the 3' end). `TrimByQuality` is Trimmomatic-style LEADING/TRAILING threshold trimming, sharing `QualityScoreAnalyzer.FindQualityTrimBounds` with `QualityScoreAnalyzer.QualityTrim`.

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

`Parse(...)` is strict in the Biopython sense (see §3.3) and accepts multi-line sequence/quality, `@`/`+` symbols inside quality strings, a repeated `+` caption, zero-length records and blank lines between records. Auto mode detects the encoding once per input, so a read confined to the overlap range (e.g. an Illumina 1.5 `BBBB` tail read) is decoded with the file's encoding. `WriteToFile` writes UTF-8 without a byte-order mark. `TrimByQuality(...)` trims only low-quality ends, not internal low-quality segments. `CalculateStatistics(...)` reports `Q20Percentage` and `Q30Percentage` as percentages, but `GcContent` as a `0..1` fraction.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- FASTQ parsing around `@` headers, `+` separator lines, and per-base quality strings.[1][2]
- Phred+33 and Phred+64 score conversion.[1][2]
- Phred/error-probability conversion formulas.[1][2]

**Intentionally simplified:**

- Auto encoding detection cannot resolve an input whose every quality character lies in ASCII 64–74 (irreducible; LIMITATIONS.md PARSE-FASTQ-001); it defaults to Phred+33 under `Permissive` and throws otherwise.
- `TrimAdapter` uses exact matching only; cutadapt's default 10 % error-tolerant alignment is not modelled.
- `Parse(TextReader)` in Auto mode buffers the records of that reader (a reader cannot be rewound); `ParseFile` and `Parse(string)` stay streaming.

**Not implemented:**

- Built-in support for FASTQ quality encodings beyond Phred+33 and Phred+64; **users should rely on:** no current alternative in this repository.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null input | Returns no records | Explicit early-return guards |
| Malformed record (no `@`, no `+`, caption mismatch, whitespace in sequence, seq/qual length mismatch, truncated) | `FormatException` | Biopython 1.88 `FastqGeneralIterator` rejects the same inputs |
| Quality symbol outside the encoding range | `FormatException` from `Parse`; `ArgumentOutOfRangeException` from `DecodeQualityScores` | Biopython `InvalidCharError` |
| Empty quality string | Detects as Phred+33 and decodes to an empty score list | Quality helpers guard empty input |
| All qualities below trim threshold | `TrimByQuality(...)` returns an empty-sequence record | End trimming can remove the full record |
| No adapter match | `TrimAdapter(...)` returns the original record | Adapter trimming is conditional |
| Adapter at position 0 | `TrimAdapter(...)` returns an empty read | cutadapt 5.2 removes the adapter and everything after it |
| Interleaved input with odd record count | Final unmatched record stays on the alternating side reached by the splitter | `SplitInterleavedReads(...)` alternates records without pair validation |

### 6.2 Limitations

Encoding detection cannot disambiguate an input confined to the Phred+33/Phred+64 overlap (ASCII 64–74), Solexa (negative-Q) FASTQ is not supported, adapter trimming is exact-match only, and `GcContent` in summary statistics is a fraction rather than a percentage.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [FastqParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/IO/FastqParserTests.cs)
- Related tests: [FastqParseTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqParseTests.cs), [FastqFilterTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqFilterTests.cs), [FastqStatisticsTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqStatisticsTests.cs), [FastqUtilityTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/FastqUtilityTests.cs)
- Test specification: [PARSE-FASTQ-001.md](../../../tests/TestSpecs/PARSE-FASTQ-001.md)

## 8. References

1. Wikipedia contributors. FASTQ format. Wikipedia. https://en.wikipedia.org/wiki/FASTQ_format
2. Cock, P.J.A., et al. 2009. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research. https://doi.org/10.1093/nar/gkp1137
3. NCBI Sequence Read Archive. FASTQ and related submit formats. https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/
4. Biopython 1.88, `Bio.SeqIO.QualityIO` (`FastqGeneralIterator`, `_get_sanger_quality_str`) — reference implementation for record validation, id splitting and write truncation.
5. FastQC, `PhredEncoding.getFastQEncodingOffset` / `PerBaseQualityScores.calculateOffsets` (https://github.com/s-andrews/FastQC) — encoding from the lowest quality character over the whole file.
6. cutadapt 5.2 (Martin 2011, EMBnet.journal 17:10) — regular 3' adapter semantics used by `TrimAdapter`.
