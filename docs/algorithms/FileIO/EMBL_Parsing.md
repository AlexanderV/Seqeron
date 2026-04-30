# EMBL Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | FileIO |
| Test Unit ID | PARSE-EMBL-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

EMBL parsing reads EMBL flat-file records used by the European Nucleotide Archive and related INSDC data distributions.[1][2] In this repository, `EmblParser` parses accession and version metadata, record descriptions, taxonomy, references, features, and sequence content from EMBL text and files. It also exposes location parsing, feature extraction, subsequence extraction, and conversion to the repository's GenBank record shape. The implementation is a core EMBL/INSDC parser rather than a full-fidelity round-trip serializer, so it preserves the main parsed fields but not every original line-level detail.[1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

EMBL flat files are line-oriented records whose lines begin with a two-character code followed by content, and whose records terminate with `//`.[1] Common line types include `ID`, `AC`, `SV`, `DE`, `KW`, `OS`, `OC`, `RN`, `RA`, `RT`, `RL`, `FT`, and `SQ`.[1] Feature locations follow INSDC syntax, including ranges, `complement(...)`, `join(...)`, `order(...)`, and partial markers using `<` and `>`.[2]

### 2.2 Core Model

Representative EMBL line types are:[1]

| Code | Description | Occurrence |
|------|-------------|------------|
| `ID` | Identification | 1 per entry |
| `AC` | Accession number(s) | One or more |
| `SV` | Sequence version | Optional |
| `DE` | Description | One or more |
| `KW` | Keywords | One or more |
| `OS` | Organism species | One or more |
| `OC` | Organism classification | One or more |
| `RN` | Reference number | One or more |
| `RA` | Reference authors | Optional |
| `RT` | Reference title | One or more |
| `RL` | Reference location | One or more |
| `FT` | Feature table | One or more |
| `SQ` | Sequence header | 1 per entry |
| `//` | Record terminator | 1 per entry |

The `ID` line uses the EMBL flat-file model documented in the original file:[1]

```text
ID   <accession>; SV <version>; <topology>; <molecule>; <class>; <division>; <length> BP.
```

Feature table entries follow the INSDC feature-table pattern:[1][2]

```text
FT   <key>           <location>
FT                   /qualifier="value"
```

Location forms preserved from the current document are:[2]

| Pattern | Description | Example |
|---------|-------------|---------|
| `n..m` | Closed range | `100..200` |
| `n` | Single position | `467` |
| `<n..m` | Partial 5' start | `<1..200` |
| `n..>m` | Partial 3' end | `100..>500` |
| `complement(loc)` | Reverse-complement strand | `complement(100..200)` |
| `join(loc1,loc2,...)` | Joined multi-part location | `join(1..50,60..100)` |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | EMBL records are terminated by `//` | EMBL flat-file entries are record-delimited by the terminator line.[1] |
| INV-02 | `complement(...)` denotes reverse-strand orientation | This is part of the INSDC location syntax shared by EMBL feature tables.[2] |
| INV-03 | `join(...)` denotes a multi-part location whose parts define an overall span from minimum start to maximum end | Joined INSDC locations are composed from ordered parts.[2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `filePath` | `string` | required | File path passed to `ParseFile(...)` | Missing or empty paths yield no records |
| `content` | `string` | required | EMBL text passed to `Parse(...)` | Null or empty input yields no records |
| `locationStr` | `string` | required | Raw location string passed to `ParseLocation(...)` | Empty input yields a zeroed location |
| `record` | `EmblRecord` | required | Parsed record passed to extraction and feature helpers | Must already contain parsed sequence and features |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Accession` | `string` | Primary accession |
| `SequenceVersion` | `string` | Sequence version from `ID` or fallback `SV` |
| `DataClass` | `string` | EMBL data class |
| `MoleculeType` | `string` | EMBL molecule type |
| `Topology` | `string` | `linear` or `circular` when present |
| `TaxonomicDivision` | `string` | EMBL taxonomic division |
| `SequenceLength` | `int` | Length from the `ID` line |
| `Description` | `string` | Joined `DE` content |
| `Keywords` | `IReadOnlyList<string>` | Parsed keyword list |
| `Organism` | `string` | Parsed organism name |
| `OrganismClassification` | `IReadOnlyList<string>` | Parsed taxonomy hierarchy |
| `References` | `IReadOnlyList<Reference>` | Parsed reference entries |
| `Features` | `IReadOnlyList<Feature>` | Parsed feature table entries |
| `Sequence` | `string` | Sequence letters extracted from the `SQ` section |
| `AdditionalFields` | `IReadOnlyDictionary<string, string>` | Non-consumed or non-standard line groups |

### 3.3 Preconditions and Validation

`ParseFile(...)` and `Parse(...)` return no records on null, empty, or missing input. Record parsing starts only on text blocks that begin with `ID` after splitting on `\n//`. When `SequenceVersion` is absent from `ID`, the parser falls back to the `SV` line; when `Accession` is absent from `ID`, it falls back to the first `AC` line. `ParseLocation("")` returns a zeroed location with empty parts. Sequence extraction keeps only letters from the `SQ` section and uppercases them.

## 4. Algorithm

### 4.1 High-Level Steps

1. Split the input into EMBL records on the `//` terminator.
2. For each record beginning with `ID`, group lines by their two-character prefixes.
3. Parse the `ID` line into accession, version, topology, molecule type, data class, division, and length.
4. Parse description, keywords, organism, classification, references, and feature-table lines.
5. Extract sequence letters from the `SQ` section.
6. Yield an `EmblRecord` and expose helper methods for location parsing, sequence extraction, and conversion.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The implementation recognizes the controlled vocabularies embedded in source for `mol_type`, EMBL data classes, and EMBL division codes while parsing `ID` fields.[1] It also preserves feature-table qualifiers as a flat key/value dictionary, with bare qualifiers represented as `"true"`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Parse` / `ParseFile` | `O(n)` | `O(n)` | Linear in record text length |
| `ParseLocation` | `O(k)` | `O(k)` | `k` = location string length |
| `ExtractSequence` | `O(m)` | `O(m)` | `m` = extracted sequence length |
| `ToGenBank` | `O(f + r)` | `O(f + r)` | `f` = feature count, `r` = reference count |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation locations:** [EmblParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/EmblParser.cs), [SequenceFormatHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/SequenceFormatHelper.cs), [FeatureLocationHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/FeatureLocationHelper.cs)

- `EmblParser.ParseFile(string)` and `EmblParser.Parse(string)`: Parse EMBL records from files or text.
- `EmblParser.ParseLocation(string)`: Parses EMBL/INSDC location syntax.
- `EmblParser.GetFeatures(...)`, `GetCDS(...)`, `GetGenes(...)`: Feature selection helpers.
- `EmblParser.ExtractSequence(...)`: Extracts subsequences for parsed feature locations.
- `EmblParser.ToGenBank(...)`: Converts an `EmblRecord` into the repository's GenBank record shape.
- `SequenceFormatHelper.ParseLocationParts(...)`: Shared GenBank/EMBL location parser.
- `FeatureLocationHelper.ExtractSequence(...)`: Shared joined/complement feature extraction helper.

### 5.2 Current Behavior

The parser splits records specifically on `\n//` and only parses trimmed blocks that begin with `ID`. `GroupLinesByPrefix(...)` concatenates repeated lines of the same prefix with spaces before field-specific parsing. `ParseFeaturesFromLines(...)` is the active feature parser and preserves EMBL `FT` indentation semantics closely enough to separate feature keys from qualifier lines. Non-consumed line groups, including standard prefixes such as `DT`, `DR`, `CC`, and `OG`, are stored in `AdditionalFields`. Feature subsequence extraction is delegated to `FeatureLocationHelper`, which handles joined locations and reverse complements and uses 1-based inclusive coordinates.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Line-prefix-based EMBL flat-file parsing with `ID`, `AC`, `SV`, `DE`, `KW`, `OS`, `OC`, `RN`, `RA`, `RT`, `RL`, `FT`, and `SQ` sections.[1]
- INSDC location parsing for ranges, `complement(...)`, `join(...)`, `order(...)`, and partial markers.[2]
- Feature-level subsequence extraction using parsed location parts.

**Intentionally simplified:**

- Record splitting depends on `\n//` delimiters and `ID`-prefixed blocks; **consequence:** malformed record separators or non-record preambles are skipped rather than repaired.
- The parser preserves feature qualifiers as a flat string dictionary; **consequence:** richer original quoting and line-layout details are not round-tripped.

**Not implemented:**

- Full EMBL occurrence-count validation for every line type and full `SQ` composition cross-checking; **users should rely on:** no current alternative in this repository.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null input | Returns no records | Explicit early-return guards |
| Minimal record with only `ID` and `SQ` | Parsed successfully | Tests cover minimal EMBL parsing |
| Missing `SV` in `ID` | Version falls back to the separate `SV` line | Repository fallback behavior |
| Empty location string | Returns a zeroed `Location` | Explicit helper behavior |
| Sequence lines containing spaces and counts | Letters are extracted and uppercased | EMBL sequence parsing strips non-letters |

### 6.2 Limitations

The parser focuses on the main EMBL flat-file fields and does not preserve the full original line structure. It is appropriate for record ingestion and feature extraction, but not for exact source-to-source round-tripping.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EmblParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/EmblParserTests.cs)
- Related tests: [EmblParseTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/EmblParseTests.cs)
- Test specification: [PARSE-EMBL-001.md](../../../tests/TestSpecs/PARSE-EMBL-001.md)
- Related algorithms: [GenBank_Parsing.md](GenBank_Parsing.md)

## 8. References

1. EMBL-EBI. EMBL User Manual. https://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
2. INSDC. Feature Table Definition. https://www.insdc.org/files/feature_table.html
3. European Nucleotide Archive. https://www.ebi.ac.uk/ena/
