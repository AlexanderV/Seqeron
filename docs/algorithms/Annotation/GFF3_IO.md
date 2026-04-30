# GFF3 I/O

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-GFF-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

GFF3 I/O in this repository is a specification-driven helper surface for reading and writing basic feature annotations. `GenomeAnnotator.ParseGff3(...)` parses flat GFF3 records into a reduced `GenomicFeature` shape, and `GenomeAnnotator.ToGff3(...)` exports `GeneAnnotation` records back to tab-delimited GFF3 lines. The implementation follows core GFF3 field rules, percent-encoding rules, and coordinate conventions, but it is intentionally narrower than the repository's fuller `GffParser` API. It should therefore be treated as a lightweight interoperability helper rather than a complete GFF3/GTF toolkit.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

GFF3 (Generic Feature Format Version 3) is a tab-delimited interchange format for genomic annotations maintained by the Sequence Ontology project. The format uses one line per feature, with comments and directives introduced by `#` and `##` respectively. The governing specification requires 1-based inclusive coordinates, defines the legal strand and phase symbols, and uses semicolon-separated `key=value` attributes in column `9`. The specification also constrains escaping in GFF3 fields to RFC 3986-style percent encoding for reserved or control characters.

### 2.2 Core Model

Each GFF3 feature line has nine ordered columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | `seqid` | Sequence or contig identifier |
| 2 | `source` | Producing program, database, or pipeline |
| 3 | `type` | Feature type, typically a Sequence Ontology term |
| 4 | `start` | 1-based inclusive start coordinate |
| 5 | `end` | 1-based inclusive end coordinate |
| 6 | `score` | Floating-point score or `.` for undefined |
| 7 | `strand` | `+`, `-`, `.`, or `?` |
| 8 | `phase` | `0`, `1`, `2` for CDS features or `.` otherwise |
| 9 | `attributes` | Semicolon-separated `key=value` pairs |

The current document also preserves the specification rule that only required characters may be percent-encoded. For column `9`, the reserved characters explicitly called out by the current source and document set are tab, newline, carriage return, percent, control characters, semicolon, equals, ampersand, and comma.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | GFF3 coordinates are 1-based and inclusive in the external text format | This is a core rule of the GFF3 specification |
| INV-02 | Attribute keys are case-sensitive in GFF3 | The specification treats `Parent` and `parent` as different keys |
| INV-03 | Phase is required for CDS features and represented as `.` for non-CDS features | This is part of the GFF3 column-8 contract documented in the spec |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `lines` | `IEnumerable<string>` | required | Raw GFF3 lines consumed by `ParseGff3(...)` | Each data line must provide all 9 tab-delimited columns to be parsed |
| `annotations` | `IEnumerable<GeneAnnotation>` | required | Gene annotations exported by `ToGff3(...)` | Expects the repository's internal `GeneAnnotation` shape |
| `seqId` | `string` | `"seq1"` | Value written to column `1` by `ToGff3(...)` | Passed through verbatim by the lightweight exporter |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ParseGff3(...)` | `IEnumerable<GenomicFeature>` | Flat feature records parsed from GFF3 lines |
| `FeatureId` | `string` | `ID` attribute when present, otherwise an auto-generated `feature_n` identifier |
| `Type` | `string` | Column `3` feature type |
| `Start` | `int` | Parsed column `4` start coordinate, preserved in 1-based inclusive GFF3 space |
| `End` | `int` | Parsed column `5` end coordinate, preserved in 1-based inclusive GFF3 space |
| `Strand` | `char` | Parsed column `7` strand symbol |
| `Score` | `double?` | Parsed column `6` score, or `null` when the field is `.` |
| `Phase` | `int?` | Parsed column `8` phase, or `null` when the field is `.` |
| `Attributes` | `IReadOnlyDictionary<string, string>` | Parsed attribute dictionary from column `9` |
| `ToGff3(...)` | `IEnumerable<string>` | Output lines beginning with `##gff-version 3` followed by serialized feature rows |

### 3.3 Preconditions and Validation

The lightweight `GenomeAnnotator` helpers expect non-null enumerables for `lines` and `annotations`; they do not add explicit null guards around those collections. `ParseGff3(...)` skips blank lines and all `#`-prefixed lines, ignores lines with fewer than 9 tab-separated fields, parses score and phase using `.` as the null sentinel, and percent-decodes attribute values with `Uri.UnescapeDataString`. `ToGff3(...)` expects the repository's internal `GeneAnnotation.Start` values in 0-based space and writes column `4` as `Start + 1` while leaving `End` unchanged. The lightweight exporter does not validate Sequence Ontology terms, `seqId` encoding, or hierarchical parent-child consistency.

## 4. Algorithm

### 4.1 High-Level Steps

1. For parsing, iterate over input lines, skip blank or `#`-prefixed lines, and split candidate data lines on tab characters.
2. If a line has all 9 fields, parse the numeric columns, decode the attributes, choose the `ID` attribute when present, and emit a `GenomicFeature` record.
3. For export, emit `##gff-version 3` as the first line.
4. Serialize each `GeneAnnotation` as a 9-column row using `seqId`, source `.`, the annotation type, 1-based start conversion, the annotation end coordinate, strand, phase, and encoded attributes.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Context | Rule |
|---------|------|
| Parse line eligibility | Fewer than 9 tab-separated fields means the line is skipped |
| Score / phase null sentinel | `.` becomes `null` on parse |
| Phase on export | `0` for `CDS`, `.` for every other exported type |
| Attribute order on export | `ID` and `product` are emitted first, then the remaining attributes except `translation` |
| GFF3 encoding on export | Encodes tab, newline, carriage return, `%`, control characters, `;`, `=`, `&`, and `,`; leaves spaces literal |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ParseGff3(...)` | `O(n)` | `O(a)` per emitted record | `n` = total input characters or lines processed; `a` = parsed attributes for the current record |
| `ToGff3(...)` | `O(n)` | `O(1)` plus output | `n` = number of annotations and emitted attribute characters |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs), [GffParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs)

- `GenomeAnnotator.ParseGff3(IEnumerable<string>)`: Lightweight parser that returns a reduced flat feature shape.
- `GenomeAnnotator.ToGff3(IEnumerable<GeneAnnotation>, string)`: Lightweight exporter from `GeneAnnotation` records.
- `GenomeAnnotator.EncodeGff3Value(string)`: Column-9 percent-encoding helper used by the exporter.
- `GffParser.Parse(...)`, `GffParser.ParseFile(...)`, `GffParser.WriteToStream(...)`, and `GffParser.BuildGeneModels(...)`: Fuller in-repository GFF/GTF surfaces.

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- `ParseGff3(...)` skips every `#`-prefixed line, so both comments and `##gff-version` / directive lines are ignored after detection.
- The parsed `GenomicFeature` shape does not preserve `seqid` or `source`; it keeps only `FeatureId`, `Type`, `Start`, `End`, `Strand`, `Score`, `Phase`, and `Attributes`.
- Missing `ID` attributes are replaced with sequential `feature_n` identifiers.
- `ToGff3(...)` always writes source `.` and exports only from `GeneAnnotation`, not from the parsed `GenomicFeature` type.
- The exporter omits the `translation` attribute and preserves spaces literally while encoding only the characters explicitly required by the GFF3 rules implemented in `EncodeGff3Value(...)`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Parsing of 9-column GFF3-like feature rows.
- Recognition of `.` as the undefined sentinel for score and phase.
- Percent decoding of encoded attribute values on parse.
- Emission of the `##gff-version 3` header, 1-based start coordinates, and phase `0` for exported CDS features.

**Intentionally simplified:**

- The lightweight parser drops `seqid` and `source` from its returned record shape; **consequence:** `GenomeAnnotator` round trips are not lossless even when the input is valid GFF3.
- The lightweight exporter accepts only `GeneAnnotation` records and always writes source `.`; **consequence:** callers cannot use it as a general-purpose GFF3 writer for arbitrary feature records.
- The exporter omits the `translation` attribute; **consequence:** large translated payloads are not preserved in lightweight GFF3 output.

**Not implemented:**

- GTF parsing and limited broader GFF-family flat-record parsing in the `GenomeAnnotator` helper surface; **users should rely on:** [GffParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs).
- Hierarchical gene-model reconstruction, multi-parent handling, and full GFF feature modeling in the lightweight helper surface; **users should rely on:** [GffParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs).
- Fully faithful parsing of classic quoted GFF2 attribute syntax; **users should rely on:** no current fully faithful GFF2-attribute parser in this repository surface.
- FASTA-section handling in the lightweight helper surface; **users should rely on:** no current FASTA-section parser in this repository surface.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `seqid` and `source` are discarded during lightweight parse | Deviation | Parsed records cannot reproduce all original columns without external state | accepted | This is specific to the reduced `GenomicFeature` record shape |
| 2 | `seqId` is written verbatim by the lightweight exporter | Assumption | Callers must supply a GFF3-compliant identifier themselves | accepted | No escaping is applied to column `1` in `ToGff3(...)` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input | Returns no parsed records | No data lines are available |
| Comment-only or directive-only input | Returns no parsed records | All `#`-prefixed lines are skipped |
| Malformed line with fewer than 9 fields | Skipped | The lightweight parser requires all 9 GFF3 columns |
| `score = "."` | `Score = null` | `.` is the undefined sentinel |
| `phase = "."` | `Phase = null` | `.` is the undefined sentinel |
| Missing `ID` attribute | Auto-generates `feature_n` | The parser fills missing feature IDs from a running counter |
| Spaces in attribute values on export | Preserved literally | The encoder follows the rule that no other printable characters should be encoded |

### 6.2 Limitations

The `GenomeAnnotator` GFF3 helpers are intentionally narrower than the repository's full parser surface. They do not preserve all columns on parse, do not write arbitrary feature records, do not validate Sequence Ontology terms, and do not expose hierarchical relationships or GTF/GFF2 behaviors. FASTA sections are also not handled by this helper surface. For richer flat-record parsing or feature-model operations, the repository's `GffParser` API is the better match for GTF and key/value-style GFF-family records, but classic quoted GFF2 attributes and FASTA-section handling are still outside the currently implemented parser family.

## 8. References

1. Stein L. Generic Feature Format Version 3 (GFF3). Sequence Ontology Project, version 1.26, 2020.
2. Wikipedia contributors. General feature format. https://en.wikipedia.org/wiki/General_feature_format
3. RFC 3986. Uniform Resource Identifier (URI): Generic Syntax.
