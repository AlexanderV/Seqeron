# TestSpec: PARSE-FASTA-001 - FASTA Parsing

## Test Unit Information
- **ID:** PARSE-FASTA-001
- **Area:** FileIO
- **Class:** FastaParser
- **Canonical Methods:** Parse, ParseFile, ParseFileAsync, ToFasta, WriteFile

## Test Classification

### Canonical vs Wrapper/Delegate
- **Canonical:** `FastaParser.Parse`, `FastaParser.ToFasta`, `FastaParser.ParseFile`, `FastaParser.WriteFile`
- **MCP Wrappers:** `ParsersTools.FastaParse`, `ParsersTools.FastaFormat`, `ParsersTools.FastaWrite`
  - Wrappers require only smoke tests (delegate to canonical)

## Test Requirements

### Must Tests (Evidence-backed)

| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M1 | Parse_SingleSequence_ReturnsCorrectEntry | Wikipedia/NCBI | Basic format parsing |
| M2 | Parse_MultiSequence_ReturnsAllEntries | Wikipedia | Multi-FASTA support |
| M3 | Parse_MultilineSequence_ConcatenatesLines | Wikipedia | Interleaved format |
| M4 | Parse_EmptyInput_ReturnsEmpty | Common edge case | Robustness |
| M5 | Parse_HeaderWithDescription_ParsesBoth | Wikipedia | Header structure |
| M6 | Parse_HeaderWithoutDescription_ParsesId | Wikipedia | Optional description |
| M7 | Parse_WhitespaceInSequence_IgnoresWhitespace | Wikipedia/NCBI | "invalid characters ignored" |
| M8 | ToFasta_SingleEntry_FormatsCorrectly | Format spec | Output formatting |
| M9 | ToFasta_LongSequence_WrapsAtLineWidth | Wikipedia | 80-char line convention |
| M10 | ToFasta_NoDescription_OmitsDescription | Format spec | Header without description |
| M11 | RoundTrip_ParseAndFormat_PreservesData | Integrity | Parse→Format→Parse identical |
| M12 | Parse_SpecialCharsInHeader_PreservesChars | NCBI identifiers | Pipes, colons in IDs |
| M13 | Parse_LeadingTrailingWhitespaceInSequence_Trimmed | NCBI | Whitespace handling |

### Should Tests

| # | Test | Rationale |
|---|------|-----------|
| S1 | ToFasta_CustomLineWidth_WrapsCorrectly | Configurable output |
| S2 | ToFasta_MultipleEntries_FormatsAll | Multi-sequence output |
| S3 | ParseFile_ValidFile_ReturnsEntries | File I/O |
| S4 | WriteFile_ValidEntries_CreatesFile | File I/O |
| S5 | Parse_BlankLinesInInput_SkipsBlankLines | Common real-world data |
| S6 | FastaEntry_Header_CombinesIdAndDescription | Property behavior |
| S7 | FastaEntry_ToString_ReturnsIdAndLength | Debugging utility |

### Could Tests

| # | Test | Rationale |
|---|------|-----------|
| C1 | ParseFileAsync_ValidFile_ReturnsEntries | Async support |
| C2 | Parse_VeryLongSequence_HandlesCorrectly | Large data |
| C3 | ToFasta_LineWidthOne_WrapsEveryChar | Extreme configuration |

## Wrapper Smoke Tests (MCP Tools)

Existing tests in `Seqeron.Mcp.Parsers.Tests`:
- `FastaParseTests` - validates MCP wrapper delegates to canonical
- `FastaFormatTests` - validates MCP wrapper delegates to canonical  
- `FastaWriteTests` - validates MCP wrapper delegates to canonical

**Status:** Adequate smoke coverage for wrappers. No changes needed.

## Consolidation Plan

### Current State
1. **Canonical tests:** `FastaParserTests.cs` - 8 tests (some gaps)
2. **MCP wrapper tests:** 3 separate files with smoke tests (adequate)

### Target State
1. **Canonical tests:** Expand `FastaParserTests.cs` to comprehensive coverage
2. **MCP wrapper tests:** Keep as-is (smoke verification only)

### Actions
- Add missing Must tests to `FastaParserTests.cs`
- Reorganize existing tests with proper regions
- Remove redundancy (none detected currently)
- Ensure all Must tests are covered

## Audit of Existing Tests

| Test | Status | Notes |
|------|--------|-------|
| Parse_SingleSequence_ParsesCorrectly | Covered | M1 ✓ |
| Parse_MultipleSequences_ParsesAll | Covered | M2 ✓ |
| Parse_NoDescription_ParsesIdOnly | Covered | M6 ✓ |
| Parse_EmptyInput_ReturnsEmpty | Covered | M4 ✓ |
| Parse_MultilineSequence_ConcatenatesLines | Covered | M3 ✓ |
| ToFasta_SingleEntry_FormatsCorrectly | Covered | M8 ✓ |
| ToFasta_LongSequence_WrapsLines | Covered | M9 ✓ |
| RoundTrip_ParseAndWrite_PreservesData | Covered | M11 ✓ |

### Missing Tests (All Closed)
- ~~M5: Parse_HeaderWithDescription_ParsesBoth~~ ✅ Covered
- ~~M7: Parse_WhitespaceInSequence_IgnoresWhitespace~~ ✅ Covered
- ~~M10: ToFasta_NoDescription_OmitsDescription~~ ✅ Covered
- ~~M12: Parse_SpecialCharsInHeader_PreservesChars~~ ✅ Covered
- ~~M13: Parse_LeadingTrailingWhitespaceInSequence_Trimmed~~ ✅ Covered
- ~~S1-S7: Should tests~~ ✅ Covered
- ~~C1-C3: Could tests~~ ✅ Covered

## Open Questions

None - FASTA format is well-documented.

## Decisions

1. **Header parsing:** First whitespace-delimited token is ID, remainder is description
2. **Empty sequences:** Not yielded (header without sequence is skipped)
3. **Default line width:** 80 characters (per historical convention)
