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
        string SequenceVersion,
        string DataClass,
        string MoleculeType,
        string Topology,
        string TaxonomicDivision,
        int SequenceLength,
        string Description,
        IReadOnlyList<string> Keywords,
        string Organism,
        IReadOnlyList<string> OrganismClassification,
        IReadOnlyList<Reference> References,
        IReadOnlyList<Feature> Features,
        string Sequence,
        IReadOnlyDictionary<string, string> AdditionalFields
    );
    
    public readonly record struct Reference(
        int Number,
        string Citation,
        string Authors,
        string Title,
        string Journal,
        string CrossReference,
        string Comment,
        string Positions,
        string Group
    );
    
    public readonly record struct Feature(
        string Key,
        Location Location,
        IReadOnlyDictionary<string, string> Qualifiers
    );
    
    public readonly record struct Location(
        int Start,
        int End,
        bool IsComplement,
        bool IsJoin,
        bool IsOrder,
        bool Is5PrimePartial,
        bool Is3PrimePartial,
        IReadOnlyList<(int Start, int End)> Parts,
        string RawLocation
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
| SHLD-24 | Parse_Reference_HasPubMed | Extract RX PUBMED cross-reference | ✅ Covered |
| SHLD-25 | Parse_Reference_HasDOI | Extract RX DOI cross-reference | ✅ Covered |
| SHLD-26 | Parse_Feature_HasQualifiers | Extract FT qualifier key-value pairs | ✅ Covered |
| SHLD-27 | Parse_Organelle_ExtractsCorrectly | Extract OG line content | ✅ Covered |
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
| COUD-01 | Parse_DatabaseCrossReference_DR | Extract DR lines to AdditionalFields | ✅ Covered |
| COUD-02 | Parse_Comments_CC | Extract CC lines to AdditionalFields | ✅ Covered |
| COUD-03 | Parse_WhitespaceOnly_ReturnsEmpty | Whitespace-only content | ✅ Covered |
| COUD-04 | ParseLocation_Order | Parse order(100..200,300..400) | ✅ Covered |
| COUD-05 | Parse_MultiLineContinuation | Handle multi-line DE, RT, etc. | ✅ Covered |
| COUD-06 | Parse_EmptyKeywords | Handle "KW   ." empty keywords | ✅ Covered |
| COUD-07 | Parse_SecondaryAccessions | Primary accession extracted when secondaries present | ✅ Covered |
| COUD-08 | Parse_DateLines_DT | Extract DT lines to AdditionalFields | ✅ Covered |
| COUD-09 | ExtractSequence_ComplementLocation | Extract reverse complement | ✅ Covered |
| COUD-10 | ExtractSequence_JoinLocation | Extract joined regions | ✅ Covered |
| COUD-11 | Parse_Organelle_OG | Extract OG line to AdditionalFields | ✅ Covered |

**Removed duplicates:**
- Old COUD-01/02 (Parse_DataClass_CON/TSA) — subsumed by SPEC-03 (AllDataClasses×10).
- Old COUD-03 (Parse_AllTaxonomicDivisions) — subsumed by SPEC-04 (AllTaxonomicDivisions×15).
- Old COUD-07 (ParseLocation_ComplementJoin) — duplicate of SHLD-33.
- Old COUD-09 (ParseLocation_RemoteReference) — documented limitation LIM-2; no test (remote entry lookups unsupported).
- Old COUD-12 (Parse_AccessionRange) — EMBL AC lines use semicolons, not ranges; replaced by COUD-07 secondary accession test.

### SPEC Tests (External Source Compliance — INSDC v11.3 & EBI Release 143)

| Test ID | Test Name | Description | Status |
|---------|-----------|-------------|--------|
| SPEC-01 | Parse_IdLine_AllInsdcMolTypes(×11) | All 11 INSDC mol_type values accepted | ✅ Covered |
| SPEC-02 | Parse_IdLine_BareDnaNotRecognisedAsMolType | Bare "DNA" rejected per vocabulary | ✅ Covered |
| SPEC-03 | Parse_IdLine_AllDataClasses(×10) | All 10 data classes including STS | ✅ Covered |
| SPEC-04 | Parse_IdLine_AllTaxonomicDivisions(×15) | All 15 division codes accepted | ✅ Covered |
| SPEC-05 | Parse_IdLine_InvalidDivisionNotRecognised | Invalid "UNK" code rejected | ✅ Covered |
| SPEC-06 | Parse_Qualifier_SlashInValue_NotTruncated | /db_xref with "/" preserved | ✅ Covered |
| SPEC-07 | Parse_Qualifier_MultipleSlashesInValue | Multiple "/" in value preserved | ✅ Covered |
| SPEC-08 | Parse_Reference_CapturesPositions | RP line stored in Reference.Positions | ✅ Covered |
| SPEC-09 | Parse_Reference_CapturesGroup | RG line stored in Reference.Group | ✅ Covered |
| SPEC-10 | Parse_EbiReferenceRecord_FullBlock | Full EBI sample reference block | ✅ Covered |
| SPEC-11 | ParseLocation_SiteBetween_ParsesRange | 123^124 site notation | ✅ Covered |

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
| Location — Order | 1 test | ✅ Complete |
| Reference — DOI/PubMed | 2 tests | ✅ Complete |
| AdditionalFields — OG/DR/CC/DT | 4 tests | ✅ Complete |
| Multi-line & Edge Cases | 3 tests | ✅ Complete |
| ExtractSequence — Complement/Join | 2 tests | ✅ Complete |
| Spec Compliance Tests | 44 tests | ✅ Complete |
| **Total** | **102 tests** | **All passing** |

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
5. **PartialLocationRecord** - Record with <, > partial indicators ✅ Added
6. **ComplexLocationRecord** - Record with complement, join, nested ✅ Added
7. **ReferenceRecord** - Record with PubMed/DOI cross-references ✅ Added
8. **OrganelleRecord** - Record with OG line (mitochondrion, plastid) ✅ Added
9. **EmptyKeywordsRecord** - Record with "KW   ." ✅ Added

## Implementation Notes

### Known Behaviors from Implementation Analysis

1. **Parse(null)** returns empty IEnumerable (not exception)
2. **Parse("")** returns empty IEnumerable
3. **Sequence** is normalized to uppercase, whitespace/numbers stripped
4. **ParseLocation** supports simple range, complement, join, order
5. **ToGenBank** converts record preserving all fields
6. **GetCDS/GetGenes** are filters on Features list

### Edge Cases Identified

1. Multi-line continuation for DE, KW, OS, OC, RA, RT, RL
2. Qualifier continuation in FT lines (col 22+)
3. Empty keyword list ("KW   .")
4. Records without sequence header summary
5. Ambiguous bases (N, R, Y, etc.) in sequence
6. Very large records (performance)

## Deviations and Assumptions

### Sources

| Source | Version | Authority |
|--------|---------|-----------|
| EBI EMBL User Manual | Release 143, March 2020 | European Nucleotide Archive — official EMBL format specification |
| INSDC Feature Table Definition | v11.3, October 2024 | International Nucleotide Sequence Database Collaboration — feature key/qualifier/location specification |

### Fixed Deviations

All deviations below were identified by systematic audit against the authoritative sources and fixed in EmblParser.cs.

| ID | Before | After | Source |
|----|--------|-------|--------|
| DEV-1 | Molecule type accepted only 5 hardcoded values (`"DNA"`, `"RNA"`, `"mRNA"`, `"genomic DNA"`, `"genomic RNA"`). Bare `"DNA"` and `"RNA"` are **not** in the INSDC controlled vocabulary. Missing 6 legitimate values. | Strict validation against the complete INSDC `/mol_type` controlled vocabulary (11 values): `"genomic DNA"`, `"genomic RNA"`, `"mRNA"`, `"tRNA"`, `"rRNA"`, `"other RNA"`, `"other DNA"`, `"transcribed RNA"`, `"viral cRNA"`, `"unassigned DNA"`, `"unassigned RNA"`. Uses `HashSet<string>` lookup. | INSDC Feature Table v11.3 — `/mol_type` qualifier |
| DEV-2 | Data class list missing `"STS"` (Sequence Tagged Site). | Added `"STS"` to `ValidDataClasses` set. All 10 data classes now covered: CON, PAT, EST, GSS, HTC, HTG, WGS, TSA, STS, STD. | EBI User Manual §3.1 |
| DEV-3 | Division validation used `trimmed.Length == 3 && char.IsUpper(trimmed[0])` — accepted any 3-char uppercase-leading string (e.g., `"UNK"`, `"BCT"`, `"ABC"`). | Strict validation against the 15 known taxonomic division codes using `HashSet<string>`: PHG, ENV, FUN, HUM, INV, MAM, VRT, MUS, PLN, PRO, ROD, SYN, TGN, UNC, VRL. | EBI User Manual §3.2 |
| DEV-4 | Qualifier parsing used regex `[^/]+` for value capture — truncated values containing `/`. E.g., `/db_xref="UniProtKB/Swiss-Prot:P26204"` yielded `"UniProtKB"` instead of full value. | Manual parsing with `IndexOf('=', 1)`: takes everything after `=` as the value, correctly preserving `/` characters in quoted values. | INSDC Feature Table v11.3 §6.2 — qualifier value delimited by `"`, not `/` |
| DEV-5 | RP (Reference Positions) line parsed as standard prefix but content discarded. | Added `Positions` field to `Reference` record struct. RP content stored in `Reference.Positions`. | EBI User Manual §3.4.10 — RP ≥1 per entry |
| DEV-6 | RG (Reference Group) line parsed as standard prefix but content discarded. | Added `Group` field to `Reference` record struct. RG content stored in `Reference.Group`. | EBI User Manual §3.4.10 — RG ≥0 per entry |
| DEV-7 | Test data used `"UNK"` and `"BCT"` as division codes. Neither exist in the EBI specification. | `"UNK"` → `"UNC"` (Unclassified). `"BCT"` → `"PRO"` (Prokaryote). | EBI User Manual §3.2 |
| DEV-8 | Test data used bare `"DNA"` as molecule type. Not in the INSDC controlled vocabulary. | `"DNA"` → `"genomic DNA"` in all test records. | INSDC Feature Table v11.3 — `/mol_type` qualifier |
| DEV-9 | OG, DR, CC, DT line content silently dropped. `IsStandardPrefix` returned true for these prefixes, excluding them from `AdditionalFields`, but no dedicated parser consumed them. | Changed AdditionalFields loop to use `consumedPrefixes` set (only prefixes actively parsed). OG, DR, CC, DT content now stored in `AdditionalFields`. | EBI User Manual §3 — all line types should be accessible |

### Remaining Known Limitations

These are parser behaviors that differ from the full specification but are documented and accepted with rationale.

| ID | Limitation | Rationale |
|----|-----------|-----------|
| LIM-1 | Site-between location `123^124` parsed as Start=123, End=124. No distinct `IsSiteBetween` flag on the Location struct. | RawLocation preserves the `^` for downstream callers that need to distinguish site from range. Practical impact is minimal as site-between describes a zero-length insertion point. |
| LIM-2 | Remote entry locations (`J00194.1:100..202`) may parse incorrectly — the version number `1` in the accession can be matched as a position by the `(\d+)` regex in `SequenceFormatHelper`. | Remote references are rare in practice and require cross-entry lookup. The `RawLocation` field preserves the original string for downstream handling. |
| LIM-3 | Deprecated single-base-from-range notation `102.110` (single period) not handled. The regex `(\d+)(?:\.\.(\d+))?` requires double-dot `..`. | Deprecated since October 2006 per INSDC spec. Illegal for new entries. |
| LIM-4 | Escaped double quotes in qualifier values (`""`) not unescaped to `"`. The string `The ""label"" qualifier` is stored as-is with `""`. | FinishQualifier trims outer quotes but does not unescape inner pairs. Callers can unescape if needed. |
| LIM-5 | `QualifierRegex` (`@"/(\w+)(?:=([^/]+))?"`) still exists as dead code. Used only by unused `ParseQualifierString`/`ParseFeatures` methods, not by the active `ParseFeaturesFromLines` code path. | Dead code retained for potential future use or removal in a cleanup pass. |

### Assumptions

None. All implementation decisions are derived from the EBI EMBL User Manual (Release 143) and the INSDC Feature Table Definition (v11.3). No internal assumptions remain.

## Verification Criteria

### Pass Criteria
- All MUST tests pass (13 tests)
- All SHOULD tests pass (34 tests)
- All COULD tests pass (11 tests)
- All SPEC tests pass (44 tests via parameterized cases)
- No regressions in existing tests

### Quality Metrics
- All public methods have at least one test
- All controlled vocabularies exhaustively tested
- All documented edge cases covered
- All assertions use exact values (no permissive Contains/GreaterThan where exact match is possible)

## Document History
- **Created**: 2025-01-28
- **Author**: Algorithm Testing Protocol
- **Last Updated**: 2026-03-12
- **Coverage Audit**: 17 weak tests strengthened (exact values, flag assertions). 12 missing tests implemented. 5 spec-level duplicates removed. Implementation fix: OG/DR/CC/DT content now stored in AdditionalFields (DEV-9). Total: 90 → 102 tests.
- **Spec Compliance Audit**: Verified against EBI EMBL User Manual Release 143 and INSDC Feature Table v11.3. 9 deviations fixed, 5 known limitations documented, 0 assumptions remaining.
