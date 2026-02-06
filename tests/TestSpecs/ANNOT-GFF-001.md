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
| **Status** | ☑ Complete |

---

## Evidence Summary

### Authoritative Sources

1. **Sequence Ontology GFF3 Specification v1.26** (Lincoln Stein, 2020)
   - Official format definition
   - URL: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

2. **Wikipedia - General feature format**
   - Overview and column descriptions
   - URL: https://en.wikipedia.org/wiki/General_feature_format

3. **RFC 3986** - URI percent-encoding rules

### Key Specification Points

- 9 tab-delimited columns
- 1-based coordinates
- "." denotes undefined values
- Percent-encoding for special characters
- Semicolon-separated key=value attributes
- Comment lines start with `#`
- Directive lines start with `##`

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
| M11 | ParseGff3 decodes URL-encoded attribute values | Percent-encoding | RFC 3986 |
| M12 | ParseGff3 skips malformed lines (< 9 fields) | Error handling | GFF3 Spec |
| M13 | ParseGff3 auto-generates ID if missing | ID attribute | Implementation |
| M14 | ToGff3 emits version header | Version directive | GFF3 Spec |
| M15 | ToGff3 converts 0-based to 1-based coordinates | Coordinate system | GFF3 Spec |
| M16 | ToGff3 URL-encodes special characters | Percent-encoding | RFC 3986 |
| M17 | ToGff3 produces valid tab-delimited output | Column format | GFF3 Spec |
| M18 | ParseGff3 handles strand values (+, -, .) | Strand field | GFF3 Spec |

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
2. **Escaping symmetry**: Decoded on parse, encoded on export
3. **Column count**: All output lines have exactly 9 tab-separated fields

---

## Audit of Existing Tests

### GenomeAnnotatorTests.cs (GFF3 region)

| Test | Classification | Action |
|------|---------------|--------|
| `ParseGff3_ValidLine_ParsesCorrectly` | Covered (M1) | Move to canonical file |
| `ParseGff3_WithScore_ParsesScore` | Covered (M3) | Move to canonical file |
| `ParseGff3_SkipsComments` | Covered (M7) | Move to canonical file |
| `ParseGff3_ParsesAttributes` | Covered (M10, M11) | Move to canonical file |
| `ToGff3_GeneratesValidOutput` | Covered (M14, M17) | Move to canonical file |
| `ToGff3_EscapesSpecialCharacters` | Covered (M16) | Move to canonical file |

### Missing Tests

- M2 (all columns), M4 (null score), M5 (phase int), M6 (null phase)
- M8 (directives), M9 (empty lines), M12 (malformed), M13 (auto ID)
- M15 (coordinate conversion), M18 (strand values)
- S3 (translation excluded), S4 (roundtrip)

---

## Consolidation Plan

1. **Create** `GenomeAnnotator_GFF3_Tests.cs` as canonical test file
2. **Move** 6 existing tests from `GenomeAnnotatorTests.cs` GFF3 regions
3. **Remove** GFF3 tests from `GenomeAnnotatorTests.cs`
4. **Add** missing Must tests (M2, M4-M6, M8-M9, M12-M13, M15, M18)
5. **Add** relevant Should tests (S3, S4)

---

## Assumptions

None. All test rationale is traceable to GFF3 specification or RFC 3986.

---

## Open Questions

None.
