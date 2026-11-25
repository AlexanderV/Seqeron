# GFF3 I/O

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-GFF-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Complete (flat feature lines; gene-model assembly lives in `GffParser`) |
| Last Reviewed | 2026-11-25 (B11 F8–F10) |

## 1. Overview

GFF3 I/O in this repository is a specification-driven helper surface for reading and writing basic feature annotations. `GenomeAnnotator.ParseGff3(...)` parses GFF3 feature lines into `GenomicFeature` records (all nine columns, multi-valued attributes as value lists), and `GenomeAnnotator.ToGff3(...)` exports `GeneAnnotation` records to spec-valid GFF3 lines. Both follow the Sequence Ontology GFF3 v1.26 column rules, percent-encoding rules and coordinate conventions, cross-checked against GenomeTools `gt gff3validator`, gffutils and BCBio-GFF ([Evidence](../../Evidence/ANNOT-GFF-001-Evidence.md)). GTF / GFF2 and gene-model assembly are in `GffParser` (PARSE-GFF-001).

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
| `seqId` | `string` | `"seq1"` | Value written to column `1` by `ToGff3(...)` | Non-empty; characters outside `[a-zA-Z0-9.:^*$@!+_?-\|]` are percent-encoded (UTF-8 bytes) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ParseGff3(...)` | `IEnumerable<GenomicFeature>` | Flat feature records parsed from GFF3 lines |
| `Seqid` / `Source` | `string` | Columns 1 / 2, percent-decoded |
| `AttributeValues` | `IReadOnlyDictionary<string, IReadOnlyList<string>>` | Column 9 tag → value list (split on the unescaped `,`, each value decoded) |
| `FeatureId` | `string` | `ID` attribute when present, otherwise an auto-generated `feature_n` identifier |
| `Type` | `string` | Column `3` feature type |
| `Start` | `int` | Parsed column `4` start coordinate, preserved in 1-based inclusive GFF3 space |
| `End` | `int` | Parsed column `5` end coordinate, preserved in 1-based inclusive GFF3 space |
| `Strand` | `char` | Parsed column `7` strand symbol |
| `Score` | `double?` | Parsed column `6` score, or `null` when the field is `.` |
| `Phase` | `int?` | Parsed column `8` phase, or `null` when the field is `.` |
| `Attributes` | `IReadOnlyDictionary<string, string>` | Column 9 tag → decoded value text (multi-values joined with `,`) |
| `ToGff3(...)` | `IEnumerable<string>` | Output lines beginning with `##gff-version 3` followed by serialized feature rows |

### 3.3 Preconditions and Validation

`ParseGff3(null)` and `ToGff3(null)` throw `ArgumentNullException`; `ToGff3` with a null/empty `seqId` throws `ArgumentException`.

`ParseGff3(...)` skips blank lines and `#`/`##` lines, drops a trailing CR, and stops at `##FASTA` or a line beginning with `>` (implied `##FASTA`). A feature line is accepted only when it satisfies the spec's column rules (the rules `gt gff3validator` enforces): exactly 9 tab-separated columns; start/end positive integers (digits only) with start ≤ end; score `.` or a floating-point number (invariant culture); strand one of `+ - . ?`; phase `.` or `0`/`1`/`2`; column 9 `.` or `;`-separated `tag=value` tokens with non-empty tag and value and no repeated tag (empty tokens from a trailing `;` are ignored). Other lines are skipped without consuming a `feature_n` number. A CDS line with phase `.` is accepted (gt accepts it).

`ToGff3(...)` throws `ArgumentException` for an annotation with an empty type, a strand outside `+ - . ?`, a negative start, an end before its start, or the interval `[0, 0)`. It does not validate Sequence Ontology terms or parent-child consistency.

## 4. Algorithm

### 4.1 High-Level Steps

1. Parsing: iterate lines; drop a trailing CR; skip blank and `#` lines; stop at `##FASTA` / `>`.
2. Split a candidate line on tabs and validate the columns (§3.3); skip the line on any violation.
3. Decode columns 1–3; split column 9 on `;`, each token at the first `=` (tag trimmed, value verbatim), each value list on the unescaped `,`, and percent-decode tags and values; emit `GenomicFeature` with `FeatureId` = `ID` or `feature_n`.
4. Export: validate every annotation, emit `##gff-version 3`, then one row per annotation: encoded seqid, encoded source (or `.`), encoded type, start (`Start + 1`, or `s` for an empty interval `[s, s)`), end, score (or `.`), strand, phase, column 9.
5. CDS phase: group CDS rows by `GeneId`, order them 5'→3' (ascending start on `+`, descending start on `−`), set the 5'-most segment to `0` and each later segment to `(3 − (Σ preceding segment lengths) mod 3) mod 3` [1].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Context | Rule |
|---------|------|
| Parse line eligibility | All column rules of §3.3; otherwise skipped |
| Score / phase null sentinel | `.` becomes `null` on parse |
| Multi-valued attributes | Parse: `Parent=a,b` → `["a","b"]`, `Note=a%2Cb` → `["a,b"]`; export: `AttributeValues` values encoded separately and joined by `,` |
| Source / score on export | the feature's real `Source` / `Score` when present; `.` when absent |
| Phase on export | `.` for non-CDS; CDS cumulative per transcript (§4.1 step 5) |
| Column 9 on export | `ID` (GeneId), `product` (Product), then `Attributes` except `translation`/`ID`/`product` and tags also in `AttributeValues`, then `AttributeValues`; empty values omitted; `.` if nothing remains |
| Encoding on export | Every column: tab, LF, CR, `%`, control chars; column 9 tags and values also `; = & ,`; seqid every character outside `[a-zA-Z0-9.:^*$@!+_?-\|]`; spaces literal elsewhere |

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
- `GenomeAnnotator.EncodeGff3Text(string, bool)` / `EncodeGff3Seqid(string)` / `DecodeGff3(string)`: private percent-encoding helpers.
- `GffParser.Parse(...)`, `GffParser.ParseFile(...)`, `GffParser.WriteToStream(...)`, and `GffParser.BuildGeneModels(...)`: Fuller in-repository GFF/GTF surfaces.

### 5.2 Current Behavior

- `ParseGff3(...)` keeps all nine columns (`Seqid`, `Source` added by B11 F8), multi-valued attributes as lists (`AttributeValues`), and skips spec-invalid lines (B11 F9).
- `ToGff3(...)` emits only spec-valid lines (validated by `gt gff3validator` on the canonical gene and a 300-annotation escaping corpus) and writes multi-valued attributes from `GeneAnnotation.AttributeValues` (B11 F10).
- Round trip: `ParseGff3(ToGff3(x))` restores seqid, source, type, ID, product, tags and value lists; Start comes back in GFF3 space (`Start + 1`).

### 5.3 Conformance to Theory / Spec

**Implemented:** the nine columns with the spec's domains; percent-encoding/decoding in every column; seqid allowed set; multi-valued column-9 attributes; `##FASTA` and implied `>` FASTA start; zero-length features; CDS phase accumulation 5'→3' per transcript [1].

**Choices (recorded, sourced):**

- CDS with phase `.` is accepted on parse: `gt gff3validator` accepts it although the spec requires a CDS phase; export always writes it.
- Malformed lines are skipped rather than aborting the parse (gt rejects the whole file).
- Numbers are parsed strictly (positive digits, invariant float); gt's C-library leniency (`+7` → 7, `-5` reset to 1, score `1,5` → 1) is not copied.
- Tags are trimmed of surrounding whitespace (as gt, gffutils and BCBio do); values keep their spaces.
- CDS segments are grouped into transcripts by `GeneId` and the 5'-most segment is phase `0` (a 5'-partial CDS starting mid-codon cannot be expressed because `GeneAnnotation` has no per-feature phase).
- The exporter omits the `translation` attribute.

**Not in this helper surface:** GTF/GFF2, hierarchical gene-model assembly, `##sequence-region` bounds checks, SO term validation — see [GffParser.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs) for GTF and gene models.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | CDS with phase `.` accepted on parse | Choice | Such a CDS has `Phase = null` | accepted | gt gff3validator accepts it |
| 2 | The 5'-most CDS segment of each transcript is phase `0` | Assumption | 5'-partial CDS not expressible | accepted | Matches the SO canonical-gene example [1] |
| 3 | Empty interval `[s, s)` written as `s..s` | Spec mapping | Parses back as `Start = End = s` (GFF3 space) | accepted | Spec zero-length feature rule |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input / only `#` lines | No records | No feature lines |
| Line not 9 columns, bad strand/phase/coordinates, repeated tag, token without `=`, empty value | Skipped | Spec column rules (gt errors) |
| `##FASTA` or `>` line | Parsing stops | Spec FASTA section |
| `Parent=a,b` / `Note=a%2Cb` | `["a","b"]` / `["a,b"]` | Spec column 9 |
| `score = "."` / `phase = "."` | `null` | Undefined sentinel |
| Missing `ID` | `feature_n` | Running counter of accepted lines |
| Export without `Source`/`Score` | `.` / `.` | Defaults |
| Export multi-segment CDS | Phase accumulates 5'→3' | [1] |
| Export invalid strand, negative start, end < start, `[0,0)`, empty type | `ArgumentException` | No GFF3 representation |
| Export seqid with space or leading `>` | Percent-encoded | Spec column 1 |
| Export empty product / attribute value | Tag omitted | gt: "has no value" |

### 6.2 Limitations

No GTF/GFF2 dialects, no gene-model assembly, no `##sequence-region` bounds check and no SO term validation in this helper surface (see `GffParser` for GTF and gene models).

## 8. References

1. Stein L. Generic Feature Format Version 3 (GFF3). Sequence Ontology Project, version 1.26, 2020.
2. Wikipedia contributors. General feature format. https://en.wikipedia.org/wiki/General_feature_format
3. RFC 3986. Uniform Resource Identifier (URI): Generic Syntax.
