# VCF Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | FileIO |
| Test Unit ID | PARSE-VCF-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

VCF parsing reads Variant Call Format records and their header metadata for genomic variation exchange.[1][2] In this repository, `VcfParser` parses VCF content from files and strings, optionally returns a structured header, classifies variants, filters records, inspects genotype/sample fields, computes summary statistics, calculates transition/transversion ratios, writes records, and exposes INFO helpers. The implementation follows the VCF text model closely for common workflows, but it remains simplified because most semantic validation is deferred to callers rather than enforced during parsing.[1]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

VCF is a tab-delimited text format with a metadata header section and a data section containing one record per variant locus.[1][2] Header lines beginning with `##` declare file format and structured metadata such as INFO, FORMAT, and FILTER fields. The `#CHROM` header line defines the fixed data columns and optional sample columns.[1]

### 2.2 Core Model

Representative header lines preserved from the current document are:[1]

```text
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=q10,Description="Quality below 10">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1
```

VCF data columns are:[1][2]

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `CHROM` | string | Chromosome identifier |
| 2 | `POS` | integer | 1-based position |
| 3 | `ID` | string | Variant identifier or `.` |
| 4 | `REF` | string | Reference allele |
| 5 | `ALT` | string | One or more alternate alleles |
| 6 | `QUAL` | float or `.` | Quality score |
| 7 | `FILTER` | string | Filter state |
| 8 | `INFO` | string | INFO annotations |
| 9 | `FORMAT` | string | Optional sample format layout |
| 10+ | `SAMPLE` | string | Optional sample values |

Genotype notation preserved from the current document is:[1]

| Notation | Meaning |
|----------|---------|
| `0/0` | Homozygous reference |
| `0/1` | Heterozygous |
| `1/1` | Homozygous alternate |
| `1/2` | Heterozygous with two alternate alleles |
| `./.` | Missing genotype |
| `0|1` | Phased heterozygous |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `POS` is 1-based | VCF coordinates are 1-based by specification.[1] |
| INV-02 | VCF data lines have at least 8 mandatory columns | The first 8 columns are fixed in the format.[1] |
| INV-03 | `PASS` and `.` represent different FILTER states | `PASS` means filters passed, while `.` means no filters were applied.[1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `filePath` | `string` | required | Path passed to `ParseFile(...)` or `ParseFileWithHeader(...)` | Missing or empty paths yield no records |
| `[Parse] content` | `string` | required | VCF text passed to `Parse(...)` | Null or empty input yields no records |
| `[ParseWithHeader] content` | `string` | required | VCF text passed to `ParseWithHeader(...)` | Null input throws; empty content yields a default header and no records |
| `records` | `IEnumerable<VcfRecord>` | required | Parsed records used by filters and statistics helpers | Enumeration may be traversed multiple times |
| `header` | `VcfHeader?` | optional | Header passed to writing helpers | Defaults to `VCFv4.3` output when omitted |
| `sampleIndex` | `int` | required | Sample index passed to genotype helpers | Out-of-range access yields null or false |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Chrom` | `string` | Chromosome identifier |
| `Pos` | `int` | 1-based variant position |
| `Id` | `string` | Parsed variant identifier |
| `Ref` | `string` | Reference allele |
| `Alt` | `string[]` | Alternate alleles |
| `Qual` | `double?` | Parsed quality score |
| `Filter` | `string[]` | Parsed filter values |
| `Info` | `IReadOnlyDictionary<string, string>` | INFO key/value pairs, with flags stored as `true` |
| `Format` | `string[]?` | Sample field layout |
| `Samples` | `IReadOnlyList<IReadOnlyDictionary<string, string>>?` | Sample dictionaries keyed by FORMAT entries |
| `VcfHeader` | `VcfHeader` | Structured `fileformat`, INFO, FORMAT, FILTER, sample names, and other metadata |

### 3.3 Preconditions and Validation

`Parse(...)` returns no records for null or empty content. `ParseWithHeader(...)` requires non-null content; empty content yields a default `VCFv4.3` header with no records, while null content throws during reader construction. `ParseFile(...)` returns no records for null, empty, or missing paths. `ParseFileWithHeader(...)` returns a default `VCFv4.3` header with empty collections and no records when the path is missing. A data line is rejected when it has fewer than 8 tab-separated columns or when `POS` is not an integer. `QUAL` becomes `null` when its value is `.`, `FILTER` becomes an empty array when its value is `.`, and INFO flags without `=` are stored as `"true"`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Read header metadata lines beginning with `##`.
2. Parse structured INFO, FORMAT, and FILTER metadata and capture other `##` lines as raw metadata.
3. Parse the `#CHROM` line to determine sample names.
4. Split each data line on tabs and parse the first 8 mandatory columns.
5. Parse optional FORMAT/sample columns when sample names are available.
6. Expose classification, filtering, genotype, statistics, Ti/Tv, writing, and INFO helper operations over the parsed records.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Variant classification preserved from the current document and aligned to repository behavior is:

| Type | Condition | Example |
|------|-----------|---------|
| `SNP` | `len(REF) = 1` and `len(ALT) = 1` | `A -> G` |
| `MNP` | `len(REF) = len(ALT) > 1` | `AT -> GC` |
| `Insertion` | `len(REF) = 1` and `len(ALT) > 1` | `A -> ATG` |
| `Deletion` | `len(REF) > 1` and `len(ALT) = 1` | `ATG -> A` |
| `Symbolic` | `ALT` is `*`, starts with `<`, or contains `[` or `]` | `A -> <DEL>` |
| `Complex` | Any other unequal-length non-symbolic case | `ATG -> CT` |

The transition/transversion ratio is computed over SNP alleles using the current document's definition:

$$
Ti/Tv = \frac{count(AG, GA, CT, TC)}{count(AC, CA, AT, TA, GC, CG, GT, TG)}
$$

Zygosity helpers use the `GT` sample field and split genotypes on `/` or `|` to distinguish homozygous reference, homozygous alternate, and heterozygous calls.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Parse` / `ParseWithHeader` | `O(n)` | `O(n)` | Linear in file text size |
| `ClassifyVariant` / `GetVariantLength` | `O(1)` | `O(1)` | Operates on one allele string |
| `FilterByChrom` / `FilterByRegion` / `FilterByQuality` / `FilterPassing` | `O(r)` | `O(m)` | `r` = number of records, `m` = matches |
| `CalculateStatistics` / `CalculateTiTvRatio` | `O(r)` | `O(1)` auxiliary | One pass over parsed records |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [VcfParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/VcfParser.cs)

- `VcfParser.ParseFile(string)`, `Parse(string)`, `ParseFileWithHeader(string)`, and `ParseWithHeader(string)`: Parse VCF records and optional header metadata.
- `VcfParser.ClassifyVariant(...)`, `IsSNP(...)`, `IsIndel(...)`, `GetVariantLength(...)`: Variant-type helpers.
- `VcfParser.FilterByChrom(...)`, `FilterByRegion(...)`, `FilterByQuality(...)`, `FilterPassing(...)`, `FilterByType(...)`, `FilterSNPs(...)`, `FilterIndels(...)`, `FilterByInfo(...)`: Filtering helpers.
- `VcfParser.GetGenotype(...)`, `IsHomRef(...)`, `IsHomAlt(...)`, `IsHet(...)`, `GetAlleleDepth(...)`, `GetReadDepth(...)`, `CalculateAlleleFrequency(...)`: Genotype and sample helpers.
- `VcfParser.CalculateStatistics(...)`, `CalculateTiTvRatio(...)`, `WriteToStream(...)`, `WriteToFile(...)`: Statistics and writing helpers.
- `VcfParser.GetInfoValue(...)`, `GetInfoInt(...)`, `GetInfoDouble(...)`, `HasInfoFlag(...)`: INFO access helpers.

### 5.2 Current Behavior

`ParseFileWithHeader(...)` returns a default `VCFv4.3` header and no records when the file path is missing. `ParseHeader(...)` structures INFO, FORMAT, and FILTER metadata and stores other `##key=value` metadata in `OtherMetadata`. Sample fields are parsed only when a `#CHROM` header defines sample names. `ClassifyVariant(...)` treats spanning deletions (`*`), symbolic alleles (`<...>`), and breakend notation containing `[` or `]` as `Symbolic`. `FilterPassing(...)` returns only records whose FILTER column is exactly `PASS`; records with `.` are excluded because they represent a different VCF state. `WriteToStream(...)` always writes a `##fileformat=` line and formats `QUAL` with two decimal places.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Parsing of VCF header metadata, `#CHROM` sample layout, and 8 mandatory data columns.[1]
- Multi-allelic ALT parsing and INFO key/value handling.[1]
- Distinct FILTER treatment for `PASS` versus `.`.[1]

**Intentionally simplified:**

- The parser stores INFO values and sample fields as strings without enforcing header-declared `Number` and `Type`; **consequence:** semantic inconsistencies remain caller-visible.
- Only INFO, FORMAT, and FILTER metadata are promoted into specialized header structures; **consequence:** other metadata remains as raw key/value strings.

**Not implemented:**

- Full allele-content validation and support for compressed BCF or BGZF-indexed workflows; **users should rely on:** no current alternative in this repository.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty content in `Parse(...)` | Returns no records | Explicit early-return guard |
| Empty content in `ParseWithHeader(...)` | Returns a default header and no records | The reader exhausts immediately without header or data lines |
| Null content in `ParseWithHeader(...)` | Throws | `StringReader` requires non-null content |
| Missing `QUAL` (`.`) | Stored as `null` | Repository normalization behavior |
| Missing `FILTER` (`.`) | Stored as an empty array | Distinct state from `PASS` |
| Data line with fewer than 8 columns | Skipped | Minimum structural validation |
| Non-integer `POS` | Skipped | Mandatory numeric coordinate parsing |
| INFO flag without `=` | Stored with value `true` | Flag-style INFO handling |

### 6.2 Limitations

The parser favors broad interoperability over strict validation. INFO and sample values remain stringly typed, semantic checks against header declarations are not enforced, and the writer normalizes output formatting rather than preserving the exact original text.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [VcfParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/VcfParserTests.cs)
- Related tests: [VcfParseTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/VcfParseTests.cs), [VcfUtilityTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/VcfUtilityTests.cs)
- Test specification: [PARSE-VCF-001.md](../../../tests/TestSpecs/PARSE-VCF-001.md)

## 8. References

1. SAMtools. HTS-specs: VCF Specification. https://samtools.github.io/hts-specs/
2. Wikipedia contributors. Variant Call Format. Wikipedia. https://en.wikipedia.org/wiki/Variant_Call_Format
