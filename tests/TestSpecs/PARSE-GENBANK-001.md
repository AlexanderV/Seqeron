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
| M3.6 | Partial 5' end (<n..m) | INSDC location syntax |
| M3.7 | Partial 3' end (n..>m) | INSDC location syntax |

#### M4: Sequence Extraction
**Source:** NCBI ORIGIN Section Definition

| Test | Description | Evidence |
|------|-------------|----------|
| M4.1 | Extract sequence from ORIGIN | NCBI field definition |
| M4.2 | Normalize to uppercase | Standard practice |
| M4.3 | Remove position numbers | NCBI format includes numbers |
| M4.4 | Remove whitespace | NCBI format includes spaces |
| M4.5 | Verify sequence characters (A,C,G,T only) | Standard bases |

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

### GenBankParserTests.cs (42 tests)
**Status:** Canonical test file, well-structured

| Category | Tests | Coverage | Assessment |
|----------|-------|----------|------------|
| Basic Parsing | 4 | Good | Complete |
| LOCUS Line | 2 | Partial | Needs division/date tests |
| Metadata | 6 | Good | Complete |
| Features | 5 | Partial | Needs qualifier tests |
| Sequence | 2 | Good | Complete |
| Location | 6 | Good | Complete |
| Utility | 5 | Good | Complete |
| Multi-record | 2 | Good | Complete |

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

1. **Canonical file:** `GenBankParserTests.cs`
2. **Additions needed:**
   - Division code parsing test
   - Date format parsing tests
   - Partial location tests (<, >)
   - More qualifier parsing tests
3. **Wrapper tests:** SequenceIOTests.cs - keep as-is (smoke tests)
4. **MCP tests:** GenBankParseTests.cs - keep separate

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

## Open Questions

None - format is well-documented by NCBI.

---

## Decisions

| Decision | Rationale |
|----------|-----------|
| Keep existing test structure | Already well-organized |
| Add missing division/date tests | Evidence shows these are documented |
| Add partial location tests | INSDC documents < and > syntax |
| Don't test all 18 divisions | Parsing is uniform, one test sufficient |
