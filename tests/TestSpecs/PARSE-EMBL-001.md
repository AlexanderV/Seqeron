# Test Specification: PARSE-EMBL-001 (EMBL Parsing)

## Test Unit Identification
- **Test Unit ID**: PARSE-EMBL-001
- **Algorithm Name**: EMBL Parsing
- **Category**: FileIO/Parsing
- **Status**: ☑ Complete
- **Evidence Document**: [PARSE-EMBL-001-Evidence.md](../Evidence/PARSE-EMBL-001-Evidence.md)

## Canonical Implementation

### Location
`src/Seqeron/Algorithms/Seqeron.Genomics.IO/EmblParser.cs`

### Public API Surface

```csharp
public static partial class EmblParser
{
    // Core Parsing Methods
    public static IEnumerable<EmblRecord> Parse(string content);
    public static IEnumerable<EmblRecord> ParseFile(string filePath);
    
    // Location Parsing
    public static Location ParseLocation(string locationStr);
    
    // Conversion
    public static GenBankParser.GenBankRecord ToGenBank(EmblRecord embl);
    
    // Utility Methods
    public static IEnumerable<Feature> GetFeatures(EmblRecord record, string featureKey);
    public static IEnumerable<Feature> GetCDS(EmblRecord record);
    public static IEnumerable<Feature> GetGenes(EmblRecord record);
    public static string ExtractSequence(EmblRecord record, Location location);
    
    // Record Types
    public readonly record struct EmblRecord(
        string Accession,
        int Version,
        string Topology,
        string MoleculeType,
        string DataClass,
        string Division,
        int SequenceLength,
        string Description,
        List<string> Keywords,
        string Organism,
        List<string> OrganismClassification,
        string Organelle,
        List<Reference> References,
        List<Feature> Features,
        string Sequence
    );
    
    public readonly record struct Reference(
        int Number,
        string Positions,
        string Authors,
        string Title,
        string Location,
        string PubMed,
        string DOI
    );
    
    public readonly record struct Feature(
        string Key,
        Location Location,
        Dictionary<string, string> Qualifiers
    );
    
    public readonly record struct Location(
        int Start,
        int End,
        bool IsComplement,
        bool IsJoin,
        bool IsPartialStart,
        bool IsPartialEnd,
        List<Location> Parts,
        string Raw
    );
}
```

## Canonical Test File
`tests/Seqeron/Seqeron.Genomics.Tests/EmblParserTests.cs`

## Test Categories

### MUST Tests (Critical - Required for Correctness)

| Test ID | Test Name | Description | Status |
|---------|-----------|-------------|--------|
| MUST-01 | Parse_ValidRecord_ReturnsOneRecord | Parse single valid EMBL record | ✅ Covered |
| MUST-02 | Parse_EmptyContent_ReturnsEmpty | Empty string returns empty collection | ✅ Covered |
| MUST-03 | Parse_NullContent_ReturnsEmpty | Null content returns empty collection | ✅ Covered |
| MUST-04 | Parse_IdLine_ExtractsAccession | Extract primary accession from ID line | ✅ Covered |
| MUST-05 | Parse_IdLine_ExtractsTopology | Extract topology (linear/circular) | ✅ Covered |
| MUST-06 | Parse_IdLine_ExtractsMoleculeType | Extract molecule type from ID line | ✅ Covered |
| MUST-07 | Parse_IdLine_ExtractsLength | Extract sequence length from ID line | ✅ Covered |
| MUST-08 | Parse_Sequence_ExtractsAndNormalizes | Extract and normalize sequence data | ✅ Covered |
| MUST-09 | Parse_Sequence_UpperCase | Sequence normalized to uppercase | ✅ Covered |
| MUST-10 | Parse_MultipleRecords_ParsesAll | Multiple records separated by // | ✅ Covered |
| MUST-11 | Parse_Features_ExtractsAll | Extract all feature table entries | ✅ Covered |
| MUST-12 | ParseLocation_SimpleRange_ParsesCorrectly | Parse location like "100..200" | ✅ Covered |
| MUST-13 | Parse_MinimalRecord_ParsesSuccessfully | Parse record with minimal required fields | ✅ Covered |

### SHOULD Tests (Important - Expected Functionality)

| Test ID | Test Name | Description | Status |
|---------|-----------|-------------|--------|
| SHLD-01 | Parse_Description_ExtractsCorrectly | Extract DE line content | ✅ Covered |
| SHLD-02 | Parse_Keywords_ParsesMultiple | Extract and parse KW line keywords | ✅ Covered |
| SHLD-03 | Parse_Organism_ExtractsCorrectly | Extract OS line organism name | ✅ Covered |
| SHLD-04 | Parse_OrganismClassification_ParsesHierarchy | Extract OC line taxonomy | ✅ Covered |
| SHLD-05 | Parse_References_ExtractsAll | Extract reference blocks | ✅ Covered |
| SHLD-06 | Parse_Reference_HasAuthors | Extract RA line authors | ✅ Covered |
| SHLD-07 | Parse_Reference_HasTitle | Extract RT line title | ✅ Covered |
| SHLD-08 | Parse_GeneFeature_HasCorrectKey | Identify gene features | ✅ Covered |
| SHLD-09 | Parse_CircularTopology_ParsesCorrectly | Handle circular topology records | ✅ Covered |
| SHLD-10 | ParseLocation_Complement_DetectsStrand | Parse complement(100..200) | ✅ Covered |
| SHLD-11 | ParseLocation_Join_ExtractsParts | Parse join(1..50,60..100) | ✅ Covered |
| SHLD-12 | Parse_Sequence_RemovesNumbersAndSpaces | Strip line numbers and whitespace | ✅ Covered |
| SHLD-13 | Parse_MultipleRecords_EachHasCorrectSequence | Each record has distinct sequence | ✅ Covered |
| SHLD-14 | ToGenBank_ConvertsSuccessfully | Convert EMBL to GenBank format | ✅ Covered |
| SHLD-15 | ToGenBank_PreservesFeatures | Features preserved in conversion | ✅ Covered |
| SHLD-16 | GetCDS_ReturnsOnlyCDSFeatures | Filter CDS features only | ✅ Covered |
| SHLD-17 | GetGenes_ReturnsOnlyGeneFeatures | Filter gene features only | ✅ Covered |
| SHLD-18 | Parse_IdLine_ExtractsSequenceVersion | Extract SV version number | ✅ Covered |
| SHLD-19 | Parse_IdLine_ExtractsDataClass | Extract data class (STD, EST, etc.) | ✅ Covered |
| SHLD-20 | Parse_IdLine_ExtractsTaxonomicDivision | Extract taxonomic division | ✅ Covered |
| SHLD-21 | ParseFile_ValidFile_ParsesSuccessfully | Parse from file path | ✅ Covered |
| SHLD-22 | ParseFile_InvalidPath_ReturnsEmpty | Handle missing file (returns empty) | ✅ Covered |
| SHLD-23 | Parse_Reference_HasJournalLocation | Extract RL line journal location | ✅ Covered |
| SHLD-24 | Parse_Reference_HasPubMed | Extract RX PUBMED cross-reference | ❌ Missing |
| SHLD-25 | Parse_Reference_HasDOI | Extract RX DOI cross-reference | ❌ Missing |
| SHLD-26 | Parse_Feature_HasQualifiers | Extract FT qualifier key-value pairs | ✅ Covered |
| SHLD-27 | Parse_Organelle_ExtractsCorrectly | Extract OG line content | ❌ Missing |
| SHLD-28 | ParseLocation_PartialStart_DetectsPartial | Parse <100..200 | ✅ Covered |
| SHLD-29 | ParseLocation_PartialEnd_DetectsPartial | Parse 100..>200 | ✅ Covered |
| SHLD-30 | ParseLocation_SingleBase_ParsesCorrectly | Parse single position "467" | ✅ Covered |
| SHLD-31 | GetFeatures_FiltersByKey | Filter features by arbitrary key | ✅ Covered |
| SHLD-32 | ExtractSequence_ReturnsSubsequence | Extract sequence for location | ✅ Covered |
| SHLD-33 | ParseLocation_ComplementJoin | Parse complement(join(...)) | ✅ Covered |
| SHLD-34 | Parse_Feature_GeneQualifier | Extract /gene qualifier | ✅ Covered |

### COULD Tests (Nice to Have - Extended Coverage)

| Test ID | Test Name | Description | Status |
|---------|-----------|-------------|--------|
| COUD-01 | Parse_DataClass_CON | Handle CON (constructed) records | ❌ Missing |
| COUD-02 | Parse_DataClass_TSA | Handle TSA records | ❌ Missing |
| COUD-03 | Parse_AllTaxonomicDivisions | Test all division codes | ❌ Missing |
| COUD-04 | Parse_DatabaseCrossReference_DR | Extract DR lines | ❌ Missing |
| COUD-05 | Parse_Comments_CC | Extract CC lines | ❌ Missing |
| COUD-06 | Parse_WhitespaceOnly_ReturnsEmpty | Whitespace-only content | ❌ Missing |
| COUD-07 | ParseLocation_ComplementJoin | Parse complement(join(...)) | ❌ Missing |
| COUD-08 | ParseLocation_Order | Parse order(100..200,300..400) | ❌ Missing |
| COUD-09 | ParseLocation_RemoteReference | Parse J00194.1:100..200 | ❌ Missing |
| COUD-10 | Parse_MultiLineContinuation | Handle multi-line DE, RT, etc. | ❌ Missing |
| COUD-11 | Parse_EmptyKeywords | Handle "KW   ." empty keywords | ❌ Missing |
| COUD-12 | Parse_AccessionRange | Handle AC ranges like X00001-X00005 | ❌ Missing |
| COUD-13 | Parse_SecondaryAccessions | Extract secondary accession numbers | ❌ Missing |
| COUD-14 | Parse_DateLines_DT | Extract creation and update dates | ❌ Missing |
| COUD-15 | ExtractSequence_ComplementLocation | Extract reverse complement | ❌ Missing |
| COUD-16 | ExtractSequence_JoinLocation | Extract joined regions | ❌ Missing |

## Test Audit Summary

### Existing Test Coverage (EmblParserTests.cs)

| Region | Tests | Status |
|--------|-------|--------|
| Basic Parsing Tests | 4 tests | ✅ Complete |
| ID Line Tests | 8 tests | ✅ Complete |
| Metadata Tests | 4 tests | ✅ Complete |
| Reference Tests | 4 tests | ✅ Complete |
| Feature Tests | 4 tests | ✅ Complete |
| Sequence Tests | 3 tests | ✅ Complete |
| Location Parsing Tests | 8 tests | ✅ Complete |
| Conversion Tests | 2 tests | ✅ Complete |
| Utility Method Tests | 4 tests | ✅ Complete |
| Multiple Records Tests | 2 tests | ✅ Complete |
| ParseFile Tests | 2 tests | ✅ Complete |
| Edge Case Tests | 2 tests | ✅ Complete |
| **Total** | **46 tests** | **All passing** |

### MCP Wrapper Tests (EmblParseTests.cs)
- 3 test fixtures for MCP delegate bindings
- These are smoke tests, not algorithm tests
- **Classification**: Supplementary (not counted toward coverage)

## Test Data Requirements

### Required Test Records

1. **SimpleEmblRecord** - Full valid record with all common fields ✅ Exists
2. **MinimalRecord** - Minimum valid record ✅ Exists
3. **CircularRecord** - Circular topology plasmid ✅ Exists
4. **MultipleRecords** - Three records in one string ✅ Exists
5. **PartialLocationRecord** - Record with <, > partial indicators ❌ Missing
6. **ComplexLocationRecord** - Record with complement, join, nested ❌ Missing
7. **ReferenceRecord** - Record with PubMed/DOI cross-references ❌ Missing
8. **OrganelleRecord** - Record with OG line (mitochondrion, plastid) ❌ Missing
9. **EmptyKeywordsRecord** - Record with "KW   ." ❌ Missing

## Implementation Notes

### Known Behaviors from Implementation Analysis

1. **Parse(null)** returns empty IEnumerable (not exception)
2. **Parse("")** returns empty IEnumerable
3. **Sequence** is normalized to uppercase, whitespace/numbers stripped
4. **ParseLocation** supports simple range, complement, join
5. **ToGenBank** converts record preserving all fields
6. **GetCDS/GetGenes** are filters on Features list

### Edge Cases Identified

1. Multi-line continuation for DE, KW, OS, OC, RA, RT, RL
2. Qualifier continuation in FT lines (col 22+)
3. Empty keyword list ("KW   .")
4. Records without sequence header summary
5. Ambiguous bases (N, R, Y, etc.) in sequence
6. Very large records (performance)

## Test Execution Plan

### Phase 1: Verify Existing Tests
1. Run existing EmblParserTests.cs
2. Confirm all 30 tests pass
3. Document any failures

### Phase 2: Add Missing SHOULD Tests
Priority order:
1. SHLD-18, SHLD-19, SHLD-20 (ID line complete coverage)
2. SHLD-21, SHLD-22 (ParseFile tests)
3. SHLD-26 (Feature qualifiers)
4. SHLD-28, SHLD-29, SHLD-30 (Partial locations)
5. SHLD-31, SHLD-32 (Utility methods)

### Phase 3: Add COULD Tests (Optional)
Based on time/priority assessment.

## Verification Criteria

### Pass Criteria
- All MUST tests pass (13 tests)
- All SHOULD tests pass (32 tests)
- No regressions in existing tests

### Quality Metrics
- Minimum 90% branch coverage for Parse method
- All public methods have at least one test
- All documented edge cases covered

## Document History
- **Created**: 2025-01-28
- **Author**: Algorithm Testing Protocol
- **Last Updated**: 2025-01-28
