# PARSE-GENBANK-001: GenBank Parsing TestSpec

## Test Unit Identification

| Field | Value |
|-------|-------|
| **ID** | PARSE-GENBANK-001 |
| **Title** | GenBank Flat File Parsing |
| **Area** | FileIO |
| **Canonical Class** | `GenBankParser` |
| **Complexity** | O(n) |

---

## Methods Under Test

| Method | Type | Description |
|--------|------|-------------|
| `Parse(string content)` | Canonical | Parse GenBank text content |
| `ParseFile(string filePath)` | File API | Parse GenBank file |
| `ParseLocation(string locationStr)` | Utility | Parse feature location string |

---

## Test Categories

### Must Tests (Evidence-Based)

#### M1: LOCUS Line Parsing
**Source:** NCBI Sample Record Documentation

| Test | Description | Evidence |
|------|-------------|----------|
| M1.1 | Extract locus name | NCBI field definition |
| M1.2 | Extract sequence length | NCBI field definition |
| M1.3 | Extract molecule type (DNA/RNA) | NCBI field definition |
| M1.4 | Extract topology (linear/circular) | NCBI field definition |
| M1.5 | Extract division code | NCBI field definition |
| M1.6 | Extract modification date | NCBI field definition |

#### M2: Metadata Extraction
**Source:** NCBI Sample Record Documentation

| Test | Description | Evidence |
|------|-------------|----------|
| M2.1 | Extract DEFINITION | NCBI field definition |
| M2.2 | Extract ACCESSION | NCBI field definition |
| M2.3 | Extract VERSION | NCBI field definition |
| M2.4 | Parse KEYWORDS (multiple values) | NCBI field definition |
| M2.5 | Handle empty KEYWORDS (period only) | NCBI: "." means no keywords |
| M2.6 | Extract ORGANISM name | NCBI field definition |
| M2.7 | Extract taxonomy lineage | NCBI field definition |

#### M3: Feature Location Parsing
**Source:** INSDC Feature Table Definition

| Test | Description | Evidence |
|------|-------------|----------|
| M3.1 | Simple range (n..m) | INSDC location syntax |
| M3.2 | Single position | INSDC location syntax |
| M3.3 | Complement location | INSDC: complement() |
| M3.4 | Join location with parts | INSDC: join() |
| M3.5 | Complement with join | INSDC: complement(join()) |
| M3.6 | Partial 5' end (<n..m) — `Is5PrimePartial` | INSDC 3.4.2.1 |
| M3.7 | Partial 3' end (n..>m) — `Is3PrimePartial` | INSDC 3.4.2.1 |
| M3.8 | Both ends partial (<n..>m) | INSDC 3.4.2.1 |
| M3.9 | Order operator — `IsOrder` | INSDC 3.4.2.2: order() |

#### M4: Sequence Extraction
**Source:** NCBI ORIGIN Section Definition

| Test | Description | Evidence |
|------|-------------|----------|
| M4.1 | Extract sequence from ORIGIN | NCBI field definition |
| M4.2 | Normalize to uppercase | Standard practice |
| M4.3 | Remove position numbers | NCBI format includes numbers |
| M4.4 | Remove whitespace | NCBI format includes spaces |
| M4.5 | Verify sequence characters (IUPAC nucleotide codes) | INSDC 7.4.1: A,C,G,T,M,R,W,S,Y,K,V,H,D,B,N |
| M4.6 | Parse sequence with IUPAC ambiguity codes | INSDC 7.4.1 |

#### M5: Edge Cases (Evidence-Based)
**Source:** NCBI Documentation, Implementation Analysis

| Test | Description | Evidence |
|------|-------------|----------|
| M5.1 | Empty content → empty result | Defensive handling |
| M5.2 | Null content → empty result | Defensive handling |
| M5.3 | Multiple records in file | NCBI: // delimiter |
| M5.4 | Record without features | Features optional |
| M5.5 | Record with minimal fields | Graceful degradation |

---

### Should Tests (Recommended)

#### S1: Feature Parsing
| Test | Description | Rationale |
|------|-------------|-----------|
| S1.1 | Extract gene features | Common feature type |
| S1.2 | Extract CDS features | Common feature type |
| S1.3 | Parse feature qualifiers | Standard qualifier format |
| S1.4 | Handle multi-value qualifiers | Qualifiers can have multiple values |

#### S2: Reference Parsing
| Test | Description | Rationale |
|------|-------------|-----------|
| S2.1 | Parse reference number | Reference identification |
| S2.2 | Extract AUTHORS | Common field |
| S2.3 | Extract TITLE | Common field |
| S2.4 | Extract JOURNAL | Common field |
| S2.5 | Extract PUBMED ID | Links to literature |

#### S3: Utility Methods
| Test | Description | Rationale |
|------|-------------|-----------|
| S3.1 | GetCDS returns only CDS | API contract |
| S3.2 | GetGenes returns only genes | API contract |
| S3.3 | GetQualifier existing returns value | API contract |
| S3.4 | GetQualifier missing returns null | API contract |
| S3.5 | ExtractSequence for location | Feature sequence extraction |

---

### Could Tests (Nice-to-Have)

| Test | Description | Rationale |
|------|-------------|-----------|
| C1 | Round-trip: parse → serialize → parse | Data preservation |
| C2 | Handle very large sequences | Performance |
| C3 | Handle malformed LOCUS lines | Error resilience |
| C4 | Async file parsing | API completeness |
| C5 | CDS translation | Utility feature |

---

## Existing Test Audit

### GenBankParserTests.cs (52 tests)
**Status:** Canonical test file, well-structured. All assertions use exact values.

| Category | Tests | Coverage | Assessment |
|----------|-------|----------|------------|
| Basic Parsing | 4 | Good | Complete |
| LOCUS Line | 2 | Good | Complete — all 6 fields (M1.1–M1.6) |
| Metadata | 7 | Good | Complete — includes empty keywords "." (M2.5) |
| Features | 3 | Good | Complete — count, locations, qualifiers |
| Sequence | 2 | Good | Complete |
| Location | 5 | Good | Complete — range, complement, join, complement+join, single |
| Utility | 5 | Good | Complete — counts verified |
| Multi-record | 2 | Good | Complete |
| Division | 3 | Good | Complete |
| Date | 2 | Good | Complete — exact date values |
| Partial Location + Order | 5 | Good | Complete — Is5PrimePartial, Is3PrimePartial, IsOrder |
| Qualifiers | 2 | Good | Complete — keys and values |
| Reference | 3 | Good | Complete — exact fields |
| Sequence Validation | 4 | Good | Complete — IUPAC codes |
| Edge Cases | 3 | Good | Complete |

### SequenceIOTests.cs (GenBank region)
**Status:** Higher-level SequenceIO wrapper tests

| Tests | Assessment |
|-------|------------|
| 5 GenBank tests | Smoke tests for SequenceIO wrapper |

**Decision:** Keep as smoke tests, don't duplicate in canonical file.

### Seqeron.Mcp.Parsers.Tests/GenBankParseTests.cs
**Status:** MCP binding tests

| Tests | Assessment |
|-------|------------|
| ~15 tests | MCP tool binding verification |

**Decision:** Keep separate - tests MCP layer, not core parser.

---

## Consolidation Plan

1. **Canonical file:** `GenBankParserTests.cs` — 52 tests, all categories covered
2. **Wrapper tests:** SequenceIOTests.cs — keep as-is (smoke tests)
3. **MCP tests:** GenBankParseTests.cs — keep separate

---

## Test Naming Convention

```
{Method}_{Scenario}_{ExpectedResult}
```

Examples:
- `Parse_ValidRecord_ReturnsOneRecord`
- `ParseLocation_ComplementJoin_DetectsBoth`
- `GetQualifier_NonExistent_ReturnsNull`

---

## Deviations and Assumptions

None. All behavior is evidence-backed:

| # | Behavior | Justification | Source |
|---|----------|---------------|--------|
| 1 | Location struct tracks `Is5PrimePartial` / `Is3PrimePartial` | `<` and `>` operators per INSDC spec | INSDC 3.4.2.1 |
| 2 | Location struct tracks `IsOrder` for `order()` operator | Distinct from `join()` per INSDC | INSDC 3.4.2.2 |
| 3 | Sequence validation accepts full IUPAC nucleotide base codes (A,C,G,T,M,R,W,S,Y,K,V,H,D,B,N) | Standard allows ambiguity codes | INSDC 7.4.1 |
| 4 | LOCUS molecule types include ss-RNA, ds-RNA, cRNA | Full set of molecule type abbreviations | NCBI GenBank format |
| 5 | Division code detection uses known NCBI division code set (20 codes + UNK) | Replaces heuristic `Length==3 && IsUpper` | NCBI GenBank divisions |
| 6 | `FeatureLocationHelper` uses `DnaSequence.GetReverseComplementString()` for complement | Supports IUPAC ambiguity codes; `DnaSequence` constructor rejects non-ACGT | INSDC 7.4.1, `SequenceExtensions.GetComplementBase` |
| 7 | Uppercase normalization of sequence data | Standard practice — GenBank ORIGIN section uses lowercase | NCBI Sample Record |
| 8 | `//` as record delimiter | Standard GenBank flat file record separator | NCBI: "Each record ends with //" |

---

## Open Questions

None — format is well-documented by NCBI and INSDC.

---

## Decisions

| Decision | Rationale |
|----------|-----------|
| Keep existing test structure | Already well-organized |
| Add empty keywords "." test (M2.5) | NCBI documents "." means no keywords |
| Add partial location tests | INSDC documents < and > syntax |
| Don't test all 18 divisions | Parsing is uniform, one test sufficient |
| Remove duplicate complement/join tests from Features | Already covered in Location region |
| Remove duplicate minimal record test | Merged into Parse_MinimalRecord_ParsesSuccessfully |
| Use exact values in all assertions | Prevents false positives from permissive Contains/ranges |
| Remove if-guard fallbacks from feature tests | Tests must fail if parsing fails, not silently pass |
