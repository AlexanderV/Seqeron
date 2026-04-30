# GenBank Parsing

| Field | Value |
|-------|-------|
| Algorithm Group | FileIO |
| Test Unit ID | PARSE-GENBANK-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

GenBank parsing reads NCBI GenBank flat-file records that combine nucleotide sequence, metadata, references, and annotated feature tables.[1][2][3] In this repository, `GenBankParser` parses record headers, references, features, and `ORIGIN` sequence data from strings and files; exposes location parsing and feature access helpers; extracts subsequences for annotated features; and translates CDS features to amino-acid strings. The implementation follows the main GenBank/INSDC record model but is simplified relative to full GenBank validation and exact source-to-source round-tripping.[1][2][3]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

GenBank flat files are organized as keyworded sections that start in column 1 and terminate with `//`.[1][2] Standard records include `LOCUS`, `DEFINITION`, `ACCESSION`, `VERSION`, `KEYWORDS`, `SOURCE`/`ORGANISM`, `REFERENCE`, `FEATURES`, and `ORIGIN` sections.[1] Feature locations use INSDC syntax with ranges, complements, joins, and partial markers.[3]

### 2.2 Core Model

The canonical record layout preserved from the current document is:[1]

```text
LOCUS       name         length bp    type    topology division date
DEFINITION  description
ACCESSION   accession_number
VERSION     accession.version
KEYWORDS    keywords or "."
SOURCE      organism_name
  ORGANISM  scientific_name
            taxonomy; lineage; here.
REFERENCE   n  (bases x to y)
  AUTHORS   author_list
  TITLE     title
  JOURNAL   journal
  PUBMED    pmid
FEATURES             Location/Qualifiers
     feature_key     location
                     /qualifier="value"
ORIGIN
        1 sequence_data_in_60_char_lines
//
```

The `LOCUS` line fields from the current document are:[1][2]

| Field | Position | Description |
|-------|----------|-------------|
| `LOCUS` | 1-5 | Keyword |
| Name | 13-28 | Locus name |
| Length | 30-40 | Sequence length plus `bp` or `aa` |
| Type | 45-47 | DNA, RNA, mRNA, and related molecule types |
| Topology | 56-63 | `linear` or `circular` when present |
| Division | 65-67 | Three-letter GenBank division code |
| Date | 69-79 | `DD-MMM-YYYY` |

GenBank division codes preserved from the current document are:[2]

| Code | Description |
|------|-------------|
| PRI | Primate |
| ROD | Rodent |
| MAM | Other mammalian |
| VRT | Other vertebrate |
| INV | Invertebrate |
| PLN | Plant, fungal, algal |
| BCT | Bacterial |
| VRL | Viral |
| PHG | Bacteriophage |
| SYN | Synthetic |
| UNA | Unannotated |
| EST | Expressed Sequence Tag |
| PAT | Patent |
| STS | Sequence Tagged Site |
| GSS | Genome Survey Sequence |
| HTG | High-Throughput Genomic |
| HTC | High-Throughput cDNA |
| ENV | Environmental |

Feature-location forms are:[3]

| Syntax | Example | Description |
|--------|---------|-------------|
| `n..m` | `100..200` | Closed range |
| `n` | `42` | Single position |
| `<n..m` | `<1..206` | Partial at the 5' end |
| `n..>m` | `500..>600` | Partial at the 3' end |
| `complement(loc)` | `complement(100..200)` | Reverse strand |
| `join(...)` | `join(1..50,60..100)` | Discontinuous region |
| `order(...)` | `order(1..50,60..100)` | Ordered but not explicitly joined |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | GenBank records terminate with `//` | Flat-file records are delimiter-terminated.[1] |
| INV-02 | GenBank feature locations are 1-based and inclusive | INSDC feature-table coordinates use inclusive endpoints.[3] |
| INV-03 | `complement(...)`, `join(...)`, and `order(...)` are first-class location operators | They are part of the INSDC location syntax.[3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `filePath` | `string` | required | Path passed to `ParseFile(...)` | Missing or empty paths yield no records |
| `content` | `string` | required | GenBank text passed to `Parse(...)` | Null or empty input yields no records |
| `locationStr` | `string` | required | Raw GenBank location string | Empty input yields a zeroed location |
| `record` | `GenBankRecord` | required | Parsed record used by helper methods | Must already contain sequence and features |
| `feature` | `Feature` | required | Parsed feature used by qualifier and extraction helpers | Location is interpreted with GenBank semantics |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Locus` | `string` | Record locus name |
| `SequenceLength` | `int` | Declared sequence length from `LOCUS` |
| `MoleculeType` | `string` | Molecule type from `LOCUS` |
| `Topology` | `string` | Record topology |
| `Division` | `string` | GenBank division code |
| `Date` | `DateTime?` | Parsed `LOCUS` date when recognized |
| `Definition` | `string` | Record definition |
| `Accession` | `string` | Primary accession |
| `Version` | `string` | Version string |
| `Keywords` | `IReadOnlyList<string>` | Parsed keyword list |
| `Organism` | `string` | Organism name |
| `Taxonomy` | `string` | Taxonomic lineage string |
| `References` | `IReadOnlyList<Reference>` | Parsed reference entries |
| `Features` | `IReadOnlyList<Feature>` | Parsed feature-table entries |
| `Sequence` | `string` | Sequence letters extracted from `ORIGIN` |
| `AdditionalFields` | `IReadOnlyDictionary<string, string>` | Remaining non-core sections |

### 3.3 Preconditions and Validation

`ParseFile(...)` and `Parse(...)` return no records for null, empty, or missing input. Records are split on `\n//` and parsed only when the resulting block begins with `LOCUS`. `KEYWORDS` containing a single `.` are normalized to an empty keyword list. `ParseLocation("")` returns a zeroed location with no parts. Sequence parsing removes digits and spaces from `ORIGIN` content and uppercases the remaining letters.

## 4. Algorithm

### 4.1 High-Level Steps

1. Split the input on the GenBank record terminator `//`.
2. For each block beginning with `LOCUS`, parse top-level sections by header keyword.
3. Parse the `LOCUS` line into locus name, length, molecule type, topology, division, and date.
4. Parse references and features from their indented substructure.
5. Extract sequence letters from the `ORIGIN` section.
6. Expose helper operations for location parsing, feature selection, subsequence extraction, qualifier lookup, and CDS translation.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Feature qualifiers preserve the standard GenBank feature-table style:

```text
     CDS             complement(join(123..456,789..900))
                     /gene="example"
                     /product="protein"
```

CDS translation follows this decision rule in the repository:

1. If the CDS already has a `/translation` qualifier, return that value.
2. Otherwise extract the nucleotide sequence for the CDS location.
3. Translate the extracted DNA with the built-in codon table.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Parse` / `ParseFile` | `O(n)` | `O(n)` | Linear in record text length |
| `ParseLocation` | `O(k)` | `O(k)` | `k` = location string length |
| `ExtractSequence` | `O(m)` | `O(m)` | `m` = extracted sequence length |
| `TranslateCDS` | `O(m)` | `O(m)` | Translation over extracted CDS length |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation locations:** [GenBankParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/GenBankParser.cs), [SequenceFormatHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/SequenceFormatHelper.cs), [FeatureLocationHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/FeatureLocationHelper.cs)

- `GenBankParser.ParseFile(string)` and `GenBankParser.Parse(string)`: Parse GenBank records from files or text.
- `GenBankParser.ParseLocation(string)`: Parses GenBank/INSDC feature locations.
- `GenBankParser.GetFeatures(...)`, `GetCDS(...)`, `GetGenes(...)`: Feature-selection helpers.
- `GenBankParser.ExtractSequence(...)`: Extracts subsequences for parsed feature locations.
- `GenBankParser.GetQualifier(...)`: Reads one qualifier from a parsed feature.
- `GenBankParser.TranslateCDS(...)`: Returns a `/translation` qualifier when present or translates extracted DNA with the built-in codon table.
- `SequenceFormatHelper.ParseLocationParts(...)` and `FeatureLocationHelper.ExtractSequence(...)`: Shared location and sequence helpers.

### 5.2 Current Behavior

The parser splits records on `\n//` and only processes trimmed blocks that begin with `LOCUS`. Continuation lines are merged into their owning section, so multi-line fields become single logical strings. `KEYWORDS` containing only `.` become an empty list. The organism name is parsed from `SOURCE`, while the taxonomy lineage is taken from the indented `ORGANISM` subsection. Sequence extraction removes digits and whitespace from `ORIGIN` lines and uppercases the remaining letters. `TranslateCDS(...)` first returns an existing `/translation` qualifier and otherwise uses a built-in standard codon table that maps unknown codons to `X`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Parsing of core GenBank sections including `LOCUS`, `DEFINITION`, `ACCESSION`, `VERSION`, `KEYWORDS`, `SOURCE`, `REFERENCE`, `FEATURES`, and `ORIGIN`.[1][2]
- INSDC location parsing for ranges, `complement(...)`, `join(...)`, `order(...)`, and partial markers.[3]
- Feature-level subsequence extraction using 1-based inclusive coordinates.

**Intentionally simplified:**

- Feature qualifiers are preserved as flat key/value strings assembled from continuation lines; **consequence:** original qualifier formatting and richer structure are not round-tripped.
- `TranslateCDS(...)` uses a built-in standard codon table when `/translation` is absent; **consequence:** translation does not vary by alternative genetic-code metadata.

**Not implemented:**

- Validation that the declared `LOCUS` sequence length matches the parsed `ORIGIN` sequence exactly; **users should rely on:** no current alternative in this repository.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null input | Returns no records | Explicit early-return guards |
| Missing file | Returns no records | `ParseFile(...)` checks for existence |
| Multiple records in one file | Parsed by splitting on `//` | GenBank record delimiter handling |
| `KEYWORDS   .` | Produces an empty keyword list | Repository normalization behavior |
| Sequence lines with counts and spacing | Letters are extracted and uppercased | `ORIGIN` normalization |

### 6.2 Limitations

The parser is suitable for ingesting and querying GenBank-like records, but it does not preserve exact input formatting or validate every declared field against the parsed content. Sequence case is normalized to uppercase.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenBankParserTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenBankParserTests.cs)
- Related tests: [GenBankParseTests.cs](../../../tests/Seqeron/Seqeron.Mcp.Parsers.Tests/GenBankParseTests.cs)
- Test specification: [PARSE-GENBANK-001.md](../../../tests/TestSpecs/PARSE-GENBANK-001.md)
- Related algorithms: [EMBL_Parsing.md](EMBL_Parsing.md)

## 8. References

1. NCBI. GenBank Sample Record. https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
2. NCBI. GenBank Overview. https://www.ncbi.nlm.nih.gov/genbank/
3. INSDC. Feature Table Definition. https://www.insdc.org/files/feature_table.html
