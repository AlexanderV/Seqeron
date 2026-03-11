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
1. **Canonical tests:** `FastaParserTests.cs` — 29 tests, comprehensive coverage
2. **MCP wrapper tests:** 3 separate files with smoke tests (adequate)
3. **Property tests:** `FastaRoundTripProperties.cs` — 5 property-based tests

### Actions Taken
- Strengthened M8, S2, S4 with exact line assertions (replaced `Contains`/`StartsWith` with `Is.EqualTo`)
- Strengthened C2 with content verification (not just length)
- Added M14: lowercase normalization test (evidence-backed)
- No duplicates found — all tests target distinct behaviors

## Audit of Existing Tests

All Must (M1-M14), Should (S1-S7), and Could (C1-C3) tests are implemented and passing in `FastaParserTests.cs` (29 tests total).
Additional tests: M4b (whitespace-only input), HeaderWithoutSequence, TabInHeader, SequenceExactlyLineWidth, RoundTrip_SpecialCharsInHeader.

## Deviations and Assumptions

None. All behavior is evidence-backed:

| # | Behavior | Justification | Source |
|---|----------|---------------|--------|
| 1 | Header parsing: first whitespace-delimited token is ID, remainder is description | Standard FASTA header structure | Wikipedia, NCBI |
| 2 | Header without sequence is not yielded | FASTA entry = header + sequence; entry without sequence is invalid | NCBI: "The line after the FASTA definition line begins the nucleotide sequence" |
| 3 | Default line width: 80 characters | Historical convention backed by all sources | Wikipedia, NCBI, LOC |
| 4 | Internal whitespace stripped from sequence lines | Spec-mandated behavior | Wikipedia: "Anything other than a valid character would be ignored (including spaces, tabulators)" |
| 5 | Blank lines skipped | Defensive handling of common real-world data | NCBI: "Blank lines are not allowed"; common parser practice |
| 6 | Lowercase mapped to uppercase | Handled by DnaSequence constructor | Wikipedia/NCBI: "lower-case letters are accepted and are mapped into upper-case" |
