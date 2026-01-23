# GFF3 I/O (GenomeAnnotator)

## Algorithm Overview

GFF3 (Generic Feature Format Version 3) I/O operations provide parsing and export of genomic feature annotations in the widely-adopted GFF3 format. This implementation provides a simplified interface in `GenomeAnnotator` for basic GFF3 read/write operations.

## Format Specification

### Sources

- **Primary**: [Sequence Ontology GFF3 Specification v1.26](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) (Lincoln Stein, 2020)
- **Secondary**: [Wikipedia - General feature format](https://en.wikipedia.org/wiki/General_feature_format)
- **Format Authority**: Sequence Ontology Project

### GFF3 Structure (9 Tab-Delimited Columns)

| Column | Name | Description |
|--------|------|-------------|
| 1 | seqid | Sequence/chromosome identifier |
| 2 | source | Algorithm or database that generated the feature |
| 3 | type | Feature type (constrained to SO terms) |
| 4 | start | 1-based start coordinate |
| 5 | end | 1-based end coordinate (start ≤ end) |
| 6 | score | Floating point score, or "." for undefined |
| 7 | strand | +, -, ., or ? |
| 8 | phase | 0, 1, 2 for CDS features; "." for others |
| 9 | attributes | Semicolon-separated key=value pairs |

### Key Specification Rules

1. **Coordinates**: 1-based, inclusive (both start and end)
2. **Escaping**: RFC 3986 percent-encoding required for:
   - Tab (%09), newline (%0A), carriage return (%0D), percent (%25)
   - In attributes: semicolon (%3B), equals (%3D), ampersand (%26), comma (%2C)
3. **Comments**: Lines starting with `#` are comments
4. **Directives**: Lines starting with `##` are meta-directives
5. **Version**: `##gff-version 3` should be first line
6. **Attributes**: Case-sensitive; ID, Name, Parent, Target, Gap, Derives_from, Note, Dbxref, Ontology_term reserved

## Implementation Details

### GenomeAnnotator.ParseGff3

```
Signature: IEnumerable<GenomicFeature> ParseGff3(IEnumerable<string> lines)
Complexity: O(n) where n = number of lines
```

**Behavior**:
- Skips empty lines and comment lines (starting with `#`)
- Skips directive lines (starting with `##`)
- Skips malformed lines (< 9 tab-separated fields)
- Parses score as double; "." → null
- Parses phase as int; "." → null
- URL-decodes attribute values using `Uri.UnescapeDataString`
- Auto-generates FeatureId if ID attribute missing

### GenomeAnnotator.ToGff3

```
Signature: IEnumerable<string> ToGff3(IEnumerable<GeneAnnotation> annotations, string seqId = "seq1")
Complexity: O(n) where n = number of annotations
```

**Behavior**:
- Emits `##gff-version 3` header
- Converts 0-based Start to 1-based (adds 1)
- URL-encodes special characters using `Uri.EscapeDataString`
- Excludes `translation` attribute (too large)
- Fixed phase output as 0 (simplified)

## Invariants

1. **Roundtrip consistency**: Parse(ToGff3(annotations)) should preserve core fields
2. **Coordinate conversion**: Internal 0-based → GFF3 1-based on export
3. **Escaping symmetry**: Escape on export, unescape on parse

## Edge Cases

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty input | Empty result | Specification |
| Comment lines only | Empty result | Specification |
| Malformed lines (< 9 fields) | Skipped | Specification |
| Score = "." | null | Specification |
| Phase = "." | null | Specification |
| URL-encoded characters | Decoded on parse | RFC 3986 |
| Special chars in attributes | Encoded on export | RFC 3986 |
| Missing ID attribute | Auto-generated | Implementation |

## Deviations from Full GFF3

This simplified implementation does NOT support:
- Hierarchical parent-child relationships
- Discontinuous features
- FASTA section (`##FASTA`)
- Multiple values per attribute
- Full validation against Sequence Ontology

For full GFF3/GTF support, see `GffParser` class (PARSE-GFF-001).

## References

1. Stein L. "Generic Feature Format Version 3 (GFF3)". Sequence Ontology Project. Version 1.26, August 2020.
2. Wikipedia contributors. "General feature format". Wikipedia, The Free Encyclopedia.
3. RFC 3986 - Uniform Resource Identifier (URI): Generic Syntax
