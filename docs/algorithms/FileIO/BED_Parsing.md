# BED Format Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | FileIO |
| Test Unit ID | PARSE-BED-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

BED parsing reads tab-delimited genomic interval records whose required core is `chrom`, `chromStart`, and `chromEnd`.[1][2] In this repository, `BedParser` parses BED text, files, and readers; supports common BED3, BED4, BED5, BED6, and BED12 layouts; and exposes interval filtering, merging, intersection, subtraction, block expansion, statistics, and writing helpers. The coordinate model is the UCSC BED convention: `chromStart` is 0-based and `chromEnd` is non-inclusive.[1] The implementation is simplified relative to the full BED ecosystem because it focuses on standard text BED records and selected BED12 consistency checks rather than full UCSC toolchain compatibility.[1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

BED (Browser Extensible Data) is a plain-text interval format used to describe genomic regions for browsers, annotations, and interval analysis.[1][2] The format requires three columns and defines additional optional columns up to BED12.[1] Coordinates are 0-based and half-open, so a feature length is `chromEnd - chromStart`.[1]

### 2.2 Core Model

Standard BED columns are:[1][2]

| Column | Field | Required | Description |
|--------|-------|----------|-------------|
| 1 | `chrom` | Yes | Chromosome or scaffold name |
| 2 | `chromStart` | Yes | 0-based start position |
| 3 | `chromEnd` | Yes | Non-inclusive end position |
| 4 | `name` | No | Feature name |
| 5 | `score` | No | Score in the range `0..1000` |
| 6 | `strand` | No | `+`, `-`, or `.` |
| 7 | `thickStart` | No | Thick display start |
| 8 | `thickEnd` | No | Thick display end |
| 9 | `itemRgb` | No | RGB display color |
| 10 | `blockCount` | No | Number of BED12 blocks |
| 11 | `blockSizes` | No | Comma-separated block sizes |
| 12 | `blockStarts` | No | Comma-separated block starts relative to `chromStart` |

The common named variants are BED3, BED4, BED5, BED6, and BED12.[1][2]

The BED length relationship is:

$$
Length = chromEnd - chromStart
$$

For BED12 records, block coordinates expand as:

$$
blockStart_i^{absolute} = chromStart + blockStarts_i
$$

Zero-length features are allowed when `chromStart == chromEnd`, which represents an insertion point under the BED coordinate convention.[1]

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Length = chromEnd - chromStart` | BED uses 0-based half-open coordinates.[1] |
| INV-02 | Valid BED coordinates satisfy `chromStart <= chromEnd` | Zero-length intervals are allowed, but negative length is invalid.[1] |
| INV-03 | For BED12, `blockCount` must match the lengths of `blockSizes` and `blockStarts` | The BED12 block model is defined by parallel block arrays.[1] |
| INV-04 | For BED12, the first block start is `0`, blocks do not overlap, and the last `blockStart + blockSize` equals `chromEnd - chromStart` | These are the standard BED12 block constraints.[1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `filePath` | `string` | required | Path passed to `ParseFile(...)` | Missing or empty paths yield no records |
| `content` | `string` | required | BED text passed to `Parse(...)` | Null or empty input yields no records |
| `reader` | `TextReader` | required | Streamed input passed to `Parse(...)` | Read line by line |
| `format` | `BedParser.BedFormat` | `Auto` | Requested parser mode | `Auto` uses the first non-header data line to set and enforce the expected field count; explicit non-`Auto` values do not currently force BED3/BED6/BED12-specific field counts during line parsing |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Chrom` | `string` | Chromosome or scaffold identifier |
| `ChromStart` | `int` | 0-based inclusive start coordinate |
| `ChromEnd` | `int` | 0-based exclusive end coordinate |
| `Name` | `string?` | Optional feature name |
| `Score` | `int?` | Optional score, clamped to `0..1000` by the implementation |
| `Strand` | `char?` | Optional strand character `+`, `-`, or `.` |
| `ThickStart` / `ThickEnd` | `int?` | Optional BED display fields |
| `ItemRgb` | `string?` | Optional BED display field |
| `BlockCount` | `int?` | Optional BED12 block count |
| `BlockSizes` / `BlockStarts` | `int[]?` | Optional BED12 block arrays |
| `Length` | `int` | Derived feature length `ChromEnd - ChromStart` |
| `HasBlocks` | `bool` | Whether the parsed record carries BED12 block information |

### 3.3 Preconditions and Validation

The parser skips empty lines, `track` lines, `browser` lines, and `#` comments.[1] In `Auto` mode, it counts fields from the first non-header data line and skips later data lines whose field count differs. A line is rejected when it has fewer than 3 fields, when `chromStart` or `chromEnd` cannot be parsed as integers, or when `chromStart > chromEnd`. BED12 records are additionally rejected when `blockCount` does not match the block-array lengths, when the first block start is not `0`, when the last block does not end at `chromEnd - chromStart`, or when blocks overlap. The parser first splits on tabs and falls back to whitespace splitting only when fewer than 3 tab fields are present.

## 4. Algorithm

### 4.1 High-Level Steps

1. Read the input line by line.
2. Skip blank lines and BED header/comment lines.
3. In `Auto` mode, determine the expected field count from the first data line.
4. Split the line into fields and parse the required BED3 coordinates.
5. Parse optional BED4-BED12 fields when present.
6. Apply BED12 consistency checks when block fields exist.
7. Yield a `BedRecord` for each valid line.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository recognizes these header line types and skips them during parsing:[1]

| Prefix | Meaning |
|--------|---------|
| `track ` | UCSC browser display settings |
| `browser ` | UCSC browser navigation settings |
| `#` | Comment line |

BED12 block-derived coordinates are expanded as:

```text
exonStart = chromStart + blockStarts[i]
exonEnd   = exonStart + blockSizes[i]
```

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Parse` / `ParseFile` | `O(n)` | `O(1)` auxiliary | Linear scan over input lines |
| `FilterByChrom` / `FilterByRegion` / `FilterByStrand` | `O(r)` | `O(1)` auxiliary | `r` = number of records |
| `MergeOverlapping` | `O(r log r)` | `O(r)` | Dominated by sorting |
| `Intersect` | `O(a * b)` worst case | `O(b)` | Grouped by chromosome, but still quadratic in the worst case |
| `Subtract` | `O(a * b)` worst case | `O(a + b)` | Per-record interval subtraction |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [BedParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/BedParser.cs)

- `BedParser.ParseFile(string, BedFormat)`: Parses BED records from a file path.
- `BedParser.Parse(string, BedFormat)` and `BedParser.Parse(TextReader, BedFormat)`: Parse BED records from text or a reader.
- `BedParser.FilterByChrom(...)`, `FilterByRegion(...)`, `FilterByStrand(...)`, `FilterByLength(...)`, `FilterByScore(...)`: Record-level filtering helpers.
- `BedParser.MergeOverlapping(...)`, `Intersect(...)`, `Subtract(...)`, `ExpandIntervals(...)`: Interval operations over parsed records.
- `BedParser.ExpandBlocks(...)`, `GetTotalBlockLength(...)`, `GetIntrons(...)`: BED12 block helpers.
- `BedParser.CalculateStatistics(...)`, `Sort(...)`, `CalculateCoverage(...)`, `WriteToStream(...)`: Statistics, ordering, coverage, and writing utilities.

### 5.2 Current Behavior

`BedParser.Parse(...)` uses `Auto` mode by default and therefore skips mixed-width BED files once the first non-header data line establishes an expected field count. Outside `Auto` mode, the parser still follows the actual field count present on each line rather than enforcing a declared BED3/BED4/BED6/BED12 layout. The parser clamps parsed `score` values into `0..1000`, accepts only `+`, `-`, and `.` as strand values, and leaves invalid strand tokens unset. `MergeOverlapping(...)` merges touching intervals because it treats `next.ChromStart <= current.ChromEnd` as mergeable. `ExpandIntervals(...)` swaps the meaning of upstream and downstream on negative-strand records. `CalculateCoverage(...)` yields depth change points rather than one output row per base position.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Required BED3 parsing with optional BED4-BED12 fields.[1][2]
- UCSC 0-based half-open coordinate semantics.[1]
- BED12 block validation for block-count agreement, first-block origin, terminal block extent, and non-overlap.[1]

**Intentionally simplified:**

- `Auto` mode enforces a single field count for all parsed data lines; **consequence:** mixed-field BED files are partially skipped instead of being interpreted per line.
- Optional BED display fields such as `thickStart`, `thickEnd`, and `itemRgb` are parsed syntactically but not fully semantically validated; **consequence:** some display-specific inconsistencies remain caller-visible.

**Not implemented:**

- bigBed conversion and full UCSC toolchain compatibility checks; **users should rely on:** external UCSC tooling for those workflows.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null input | Returns no records | Explicit early-return guards |
| `track`, `browser`, or comment lines | Skipped | BED header/comment handling |
| `chromStart == chromEnd` | Parsed as a zero-length feature | BED permits insertion-point intervals.[1] |
| Invalid numeric coordinates | Line skipped | Parsing requires integer `chromStart` and `chromEnd` |
| Invalid BED12 block structure | Line skipped | Implementation enforces selected BED12 constraints |

### 6.2 Limitations

The parser targets standard BED text records and selected interval utilities rather than full UCSC ecosystem behavior. Auto-detection assumes a uniform data-line width across the file, and optional BED display semantics are not exhaustively validated.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [BedParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/IO/BedParserTests.cs)
- Related tests: [BedFilterTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/BedFilterTests.cs), [BedIntersectTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/BedIntersectTests.cs), [BedMergeTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/BedMergeTests.cs)
- Test specification: [PARSE-BED-001.md](../../../tests/TestSpecs/PARSE-BED-001.md)

## 8. References

1. UCSC Genome Browser. BED format FAQ. https://genome.ucsc.edu/FAQ/FAQformat.html#format1
2. Wikipedia contributors. BED (file format). Wikipedia. https://en.wikipedia.org/wiki/BED_(file_format)
