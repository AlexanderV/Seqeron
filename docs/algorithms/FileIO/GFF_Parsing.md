# GFF/GTF Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | FileIO |
| Test Unit ID | PARSE-GFF-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

GFF parsing reads tab-delimited genome annotation records in GFF3, GTF, and legacy GFF2-style layouts.[1][2][3][4] In this repository, `GffParser` parses annotation rows from strings, files, and readers; filters records by type, sequence, and region; builds simplified gene models; writes records back out; extracts reference subsequences; and merges overlapping features. The implementation supports all three declared formats, but its `Auto` mode is intentionally conservative: it only uses the `##gff-version` directive to distinguish GFF3 from legacy GFF and otherwise defaults to GFF3.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

GFF and GTF are nine-column annotation formats used to describe genes, transcripts, exons, CDS features, and related genome annotations.[1][2][4] Coordinates are 1-based and fully closed, so a feature spans both `start` and `end`.[1][2] Column 9 stores attributes, but the encoding differs between GFF3 (`key=value`) and GTF (`key "value"`).[1][4]

### 2.2 Core Model

The shared nine-column row layout is:[1][2][4]

```text
seqid  source  type  start  end  score  strand  phase  attributes
```

Column definitions preserved from the current document are:

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `seqid` | string | Chromosome or scaffold identifier |
| 2 | `source` | string | Generating program or database |
| 3 | `type` | string | Feature type such as `gene`, `mRNA`, `exon`, or `CDS` |
| 4 | `start` | int | 1-based start position |
| 5 | `end` | int | 1-based end position, inclusive |
| 6 | `score` | float? | Confidence score, or `.` when undefined |
| 7 | `strand` | char | `+`, `-`, `.`, or `?` |
| 8 | `phase` | int? | `0`, `1`, `2`, or `.` |
| 9 | `attributes` | dict | Key-value attributes |

Attribute encodings preserved from the current document are:

| Format | Syntax | Example |
|--------|--------|---------|
| GFF3 | `key=value;key=value` | `ID=gene00001;Name=EDEN;Parent=transcript001` |
| GTF | `key "value"; key "value";` | `gene_id "ENSG00001"; gene_name "TestGene";` |

Reserved-character encodings from the current document are:[1]

| Character | Encoding |
|-----------|----------|
| Tab | `%09` |
| Newline | `%0A` |
| Carriage Return | `%0D` |
| `%` | `%25` |
| `;` | `%3B` |
| `=` | `%3D` |
| `&` | `%26` |
| `,` | `%2C` |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | GFF/GTF coordinates are 1-based and fully closed | Format specifications define inclusive coordinates.[1][2][4] |
| INV-02 | `start <= end` for valid features | The formats represent closed intervals.[1][2] |
| INV-03 | `phase` is meaningful for CDS and is one of `0`, `1`, `2`, or `.` | Phase semantics are part of the column-8 definition.[1][4] |
| INV-04 | Reserved characters in GFF3 attributes are percent-encoded | GFF3 attribute escaping follows the specification's URL-escape model.[1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `filePath` | `string` | required | Path passed to `ParseFile(...)` | Missing or empty paths yield no records |
| `content` | `string` | required | GFF or GTF text passed to `Parse(...)` | Null or empty input yields no records |
| `reader` | `TextReader` | required | Reader passed to `Parse(...)` | Parsed line by line |
| `format` | `GffParser.GffFormat` | `Auto` | Selected annotation dialect | `Auto` only uses `##gff-version` and otherwise falls back to GFF3 |
| `records` | `IEnumerable<GffRecord>` | required | Parsed records used by filter, model, and merge helpers | Enumeration may be materialized internally |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Seqid` | `string` | Sequence or chromosome identifier |
| `Source` | `string` | Source column |
| `Type` | `string` | Feature type |
| `Start` / `End` | `int` | 1-based inclusive coordinates |
| `Score` | `double?` | Optional numeric score |
| `Strand` | `char` | Strand symbol |
| `Phase` | `int?` | Optional CDS phase |
| `Attributes` | `IReadOnlyDictionary<string, string>` | Parsed attributes |
| `GeneModel` | `GeneModel` | Aggregated gene, transcript, exon, CDS, and UTR collections |

### 3.3 Preconditions and Validation

The parser skips blank lines, `##` directives, and `#` comment lines. A data line is rejected when it has fewer than 8 tab-separated fields, or when `start` or `end` cannot be parsed as integers. A `.` score becomes `null`, and a `.` phase becomes `null`. Attributes are parsed only when column 9 is present; otherwise the record receives an empty attribute dictionary.

## 4. Algorithm

### 4.1 High-Level Steps

1. Read the input line by line.
2. Skip blank lines and comments.
3. In `Auto` mode, inspect `##gff-version` directives when present.
4. Split each data row on tabs and parse columns 1 through 8.
5. Parse attributes according to the resolved format.
6. Yield `GffRecord` values and expose filtering, gene-model, writing, extraction, and merge helpers.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The repository uses these attribute-parsing rules:

| Format | Parsing rule |
|--------|--------------|
| GFF3 | Split on `;`, then split each part on the first `=`, and URL-unescape keys and values |
| GTF | Split on `;`, trim each part, split on the first space, and strip surrounding quotes from the value |
| GFF2 | Reuses the non-GTF branch in the implementation |

Gene-model construction follows the hierarchy encoded in `ID` and `Parent` attributes. The implementation groups top-level `gene` records, collects direct transcript-like children of types `mRNA`, `transcript`, or `ncRNA`, and then collects transcript children whose types are `exon`, `CDS`, or contain `utr`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Parse` / `ParseFile` | `O(n)` | `O(n)` | Linear scan over annotation rows |
| `FilterByType` / `FilterByRegion` | `O(n)` | `O(m)` | `m` = number of matches |
| `BuildGeneModels` | `O(n)` | `O(n)` | Internal indexing by `ID` and `Parent` |
| `MergeOverlapping` | `O(n log n)` | `O(n)` | Dominated by sorting |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GffParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs)

- `GffParser.ParseFile(string, GffFormat)`, `Parse(string, GffFormat)`, and `Parse(TextReader, GffFormat)`: Parse GFF, GTF, or GFF2-like records.
- `GffParser.FilterByType(...)`, `FilterBySeqid(...)`, `FilterByRegion(...)`, `GetGenes(...)`, `GetExons(...)`, `GetCDS(...)`: Filtering helpers.
- `GffParser.BuildGeneModels(...)`, `GetAttribute(...)`, `GetGeneName(...)`: Hierarchy and attribute helpers.
- `GffParser.CalculateStatistics(...)`, `WriteToStream(...)`, `WriteToFile(...)`: Statistics and output helpers.
- `GffParser.ExtractSequence(...)`, `MergeOverlapping(...)`: Sequence extraction and interval merge helpers.

### 5.2 Current Behavior

`Auto` format detection only inspects the `##gff-version` directive. If that directive is absent, the implementation still parses rows as GFF3; it does not infer GTF from the attribute syntax. `ParseLine(...)` accepts eight-column rows and attaches an empty attribute dictionary when column 9 is missing. `BuildGeneModels(...)` recognizes `gene`, `mRNA`, `transcript`, and `ncRNA` as transcript-bearing hierarchy anchors and treats any feature type containing `utr` as a UTR. When a child record lists multiple `Parent` IDs, the same child record is attached under each referenced parent. `WriteToStream(...)` emits `##gff-version 3` only for GFF3 output and formats GTF attributes with trailing semicolons.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Nine-column GFF/GTF row parsing with 1-based inclusive coordinates.[1][2][4]
- GFF3 attribute escaping and GTF quoted-attribute parsing.[1][4]
- Feature-level sequence extraction using the parsed interval and strand.

**Intentionally simplified:**

- `Auto` detection only consults `##gff-version`; **consequence:** GTF inputs must be selected explicitly when they lack that directive.
- Gene-model construction recognizes a limited hierarchy of transcript, exon, CDS, and UTR feature types; **consequence:** richer ontologies remain available only as raw `GffRecord` rows.

**Not implemented:**

- Preservation of original directives, comments, and full metadata during round-trip writing; **users should rely on:** no current alternative in this repository.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null input | Returns no records | Explicit early-return guards |
| Directive or comment lines | Skipped | Parser ignores non-data lines |
| Row with fewer than 8 columns | Skipped | Minimum structural validation |
| `.` score or phase | Stored as `null` | Repository normalization behavior |
| Multiple `Parent` IDs | Child record is attached to each parent | Hierarchy builder splits comma-separated parents |

### 6.2 Limitations

The parser focuses on core row parsing and practical annotation utilities. Auto-detection is conservative, gene-model construction covers a restricted feature hierarchy, and the writer does not preserve comments or directives beyond the generated GFF3 header.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GffParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GffParserTests.cs)
- Related tests: [GffParseTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/GffParseTests.cs), [GenomeAnnotator_GFF3_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_GFF3_Tests.cs)
- Test specification: [PARSE-GFF-001.md](../../../tests/TestSpecs/PARSE-GFF-001.md)

## 8. References

1. Sequence Ontology Project. GFF3 Specification v1.26. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
2. UCSC Genome Browser. GFF format FAQ. https://genome.ucsc.edu/FAQ/FAQformat.html
3. Wikipedia contributors. General Feature Format. Wikipedia. https://en.wikipedia.org/wiki/General_feature_format
4. Washington University. GTF2.2 Specification. http://mblab.wustl.edu/GTF22.html
