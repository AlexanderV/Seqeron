# Test Specification: ANNOT-GFF-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | ANNOT-GFF-001 |
| **Area** | Annotation |
| **Algorithm** | GFF3 I/O |
| **Canonical Methods** | `GenomeAnnotator.ParseGff3(IEnumerable<string>)`, `GenomeAnnotator.ToGff3(IEnumerable<GeneAnnotation>, string)` |
| **Complexity** | O(n) |
| **Evidence Sources** | Sequence Ontology GFF3 Spec v1.26, Wikipedia (GFF), RFC 3986 |
| **Date** | 2026-03-05 |
| **Status** | Final |

---

## Evidence Summary

### Authoritative Sources

| Source | URL/Reference | Key Information |
|--------|---------------|-----------------|
| Sequence Ontology GFF3 Specification v1.26 | https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md | 9 tab-delimited columns; 1-based coordinates; encoding rules; phase; strand |
| Wikipedia: General feature format | https://en.wikipedia.org/wiki/General_feature_format | Overview: columns, phase, meta directives |
| RFC 3986 | https://tools.ietf.org/html/rfc3986 | Percent-encoding base standard |

### Key Specification Points

- **9 tab-delimited columns** (seqid, source, type, start, end, score, strand, phase, attributes)
- **1-based coordinates**, inclusive (both start and end)
- **"."** denotes undefined values (score, phase, strand)
- **Strand**: `+`, `-`, `.`, `?` (4 values, per GFF3 Spec column 7)
- **Phase**: `0`, `1`, `2` for CDS features; `"."` for all others (NOTE 4: "The phase is REQUIRED for all CDS features")
- **Comments**: lines starting with `#`; **Directives**: lines starting with `##`
- **Attributes**: semicolon-separated `key=value` pairs in column 9; attribute names are case-sensitive

---

## Encoding Mechanism

Per GFF3 Spec v1.26, Section "Description of the Format":

> GFF3 files are nine-column, tab-delimited, plain text files. Literal use of tab, newline, carriage return, the percent (%) sign, and control characters must be encoded using RFC 3986 Percent-Encoding; **no other characters may be encoded.**

> In addition, the following characters have reserved meanings in column 9 and must be escaped when used in other contexts: `;` `=` `&` `,`

> Note that unescaped spaces are allowed within fields.

**Encoded characters (column 9 attribute values)**:

| Character | Encoding | Reason |
|-----------|----------|--------|
| Tab | %09 | General |
| Newline | %0A | General |
| CR | %0D | General |
| % | %25 | General |
| Control chars | %00–%1F, %7F | General |
| ; | %3B | Column 9 reserved |
| = | %3D | Column 9 reserved |
| & | %26 | Column 9 reserved |
| , | %2C | Column 9 reserved |
| Space | **NOT encoded** | "unescaped spaces are allowed" |

---

## Test Classification

### Must Tests (Evidence-Based)

| ID | Test | Rationale | Source |
|----|------|-----------|--------|
| M1 | ParseGff3 parses valid line correctly | Core functionality | GFF3 Spec |
| M2 | ParseGff3 extracts all 9 columns | Column structure | GFF3 Spec |
| M3 | ParseGff3 parses numeric score | Score handling | GFF3 Spec |
| M4 | ParseGff3 handles "." as null score | Undefined value | GFF3 Spec |
| M5 | ParseGff3 parses phase as integer | CDS phase | GFF3 Spec |
| M6 | ParseGff3 handles "." as null phase | Undefined phase | GFF3 Spec |
| M7 | ParseGff3 skips comment lines (`#`) | Comment handling | GFF3 Spec |
| M8 | ParseGff3 skips directive lines (`##`) | Directive handling | GFF3 Spec |
| M9 | ParseGff3 skips empty lines | Empty line handling | GFF3 Spec |
| M10 | ParseGff3 parses semicolon-separated attributes | Attribute format | GFF3 Spec |
| M11 | ParseGff3 decodes percent-encoded attribute values | Percent-encoding | RFC 3986 / GFF3 Spec |
| M12 | ParseGff3 skips malformed lines (< 9 fields) | Error handling | GFF3 Spec |
| M13 | ParseGff3 auto-generates ID if missing | ID attribute | GFF3 Spec: "ID optional for features without children" |
| M14 | ToGff3 emits version header | Version directive | GFF3 Spec |
| M15 | ToGff3 converts 0-based to 1-based coordinates | Coordinate system | GFF3 Spec |
| M16 | ToGff3 encodes only GFF3-required special characters | Encoding rules | GFF3 Spec: "no other characters may be encoded" |
| M17 | ToGff3 produces valid tab-delimited output | Column format | GFF3 Spec |
| M18 | ParseGff3 handles strand values (+, -, ., ?) | Strand field | GFF3 Spec |
| M19 | ToGff3 encodes all required column 9 characters | ;, =, &, , | GFF3 Spec |
| M20 | ToGff3 outputs phase "." for non-CDS features | Phase rules | GFF3 Spec NOTE 4 |
| M21 | ToGff3 outputs phase "0" for CDS features | Phase rules | GFF3 Spec NOTE 4 |

### Should Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | ParseGff3 handles multiple attributes per line | Common use case |
| S2 | ParseGff3 handles large file simulation | Performance |
| S3 | ToGff3 excludes translation attribute | Size concern |
| S4 | Roundtrip preserves core data | Data integrity |

### Could Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | ParseGff3 handles mixed case input | Robustness |
| C2 | ToGff3 handles empty annotations | Edge case |

---

## Invariants

1. **Coordinate conversion**: Export adds 1 to 0-based internal Start
2. **Escaping symmetry**: Decoded on parse (`Uri.UnescapeDataString`), encoded on export (`EncodeGff3Value`)
3. **Column count**: All output lines have exactly 9 tab-separated fields
4. **Phase**: CDS → `0`; non-CDS → `.`

---

## Deviations from Literature

| Aspect | Literature | Implementation | Justification |
|--------|------------|----------------|---------------|
| Phase granularity | 0, 1, or 2 for CDS | Always 0 for CDS | `GeneAnnotation` record has no Phase field; phase defaults to 0 |
| Hierarchical relationships | Parent-child feature grouping | Not supported | Simplified GFF3 I/O; full hierarchy available in `GffParser` (PARSE-GFF-001) |
| Multi-value attributes | Comma-separated values for Parent, Alias, Note, etc. | Single value only | Simplified implementation |
| FASTA section | `##FASTA` directive | Not handled | Out of scope for annotation export |
| Seqid encoding | Characters not in `[a-zA-Z0-9.:^*$@!+_?-\|]` must be escaped | seqId written raw | `seqId` parameter is application-controlled |

---

## Test File Structure

```
GenomeAnnotator_GFF3_Tests.cs
├── ParseGff3 - Basic Parsing
│   ├── M1: Valid line
│   ├── M2: All 9 columns
│   └── M18: Strand values (+, -, ., ?)
├── ParseGff3 - Score Handling
│   ├── M3: Numeric score
│   └── M4: Null score
├── ParseGff3 - Phase Handling
│   ├── M5: Phase values (0, 1, 2)
│   └── M6: Null phase
├── ParseGff3 - Comment and Directive Handling
│   ├── M7: Comments
│   ├── M8: Directives
│   └── M9: Empty lines
├── ParseGff3 - Attribute Parsing
│   ├── M10: Semicolon-separated attributes
│   ├── M11: Percent-encoded values
│   └── S1: Multiple attributes
├── ParseGff3 - Error Handling
│   ├── M12: Malformed lines
│   └── M13: Auto-generated ID
├── ToGff3 - Basic Export
│   ├── M14: Version header
│   ├── M15: 1-based coordinates
│   └── M17: Tab-delimited 9-column output
├── ToGff3 - Encoding
│   ├── M16: GFF3-compliant encoding (spaces NOT encoded)
│   ├── M19: All required characters encoded
│   ├── M20: Phase "." for non-CDS
│   ├── M21: Phase "0" for CDS
│   └── S3: Translation excluded
├── ToGff3 - Edge Cases
│   └── C2: Empty annotations
├── Roundtrip Tests
│   └── S4: Core data preserved
└── Multiple Features
    └── Order preserved
```
