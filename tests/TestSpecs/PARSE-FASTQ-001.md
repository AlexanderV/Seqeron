# TestSpec: PARSE-FASTQ-001

## Test Unit Information
- **ID:** PARSE-FASTQ-001
- **Area:** FileIO
- **Class Under Test:** `FastqParser`
- **Complexity:** O(n)

---

## Methods Under Test

| Method | Type | Deep Test | Smoke Test |
|--------|------|-----------|------------|
| `Parse(content, encoding)` | Canonical | ✓ | |
| `ParseFile(filePath, encoding)` | File I/O | ✓ | |
| `DetectEncoding(qualityString)` | Helper | ✓ | |
| `DecodeQualityScores(qualityString, encoding)` | Helper | ✓ | |
| `EncodeQualityScores(scores, encoding)` | Helper | ✓ | |
| `PhredToErrorProbability(phred)` | Utility | ✓ | |
| `ErrorProbabilityToPhred(probability)` | Utility | ✓ | |
| `FilterByQuality(records, minQ)` | Filter | ✓ | |
| `FilterByLength(records, min, max)` | Filter | ✓ | |
| `TrimByQuality(record, minQ)` | Trim | ✓ | |
| `TrimAdapter(record, adapter, minOverlap)` | Trim | ✓ | |
| `CalculateStatistics(records)` | Stats | ✓ | |
| `InterleavePairedReads(r1, r2)` | Paired | ✓ | |
| `SplitInterleavedReads(interleaved)` | Paired | ✓ | |
| `WriteToFile(filePath, records)` | File I/O | ✓ | |
| `ToFastqString(record)` | Utility | ✓ | |

---

## Test Categories

### MUST Tests (Evidence-Based)

#### M1. Basic Parsing
- **M1.1** Parse single record extracts Id, Description, Sequence, QualityString correctly
- **M1.2** Parse multiple records returns correct count
- **M1.3** Sequence length equals QualityString length for each record
- **M1.4** Empty content returns empty enumerable (no exception)
- **M1.5** Null content returns empty enumerable (no exception)

*Evidence: Wikipedia FASTQ format specification*

#### M2. Quality Encoding Detection
- **M2.1** DetectEncoding returns Phred33 when quality contains chars < '@' (ASCII 64)
- **M2.2** DetectEncoding returns Phred64 when quality contains chars > 'I' (ASCII 73)
- **M2.3** DetectEncoding defaults to Phred33 for ambiguous range
- **M2.4** DetectEncoding returns Phred33 for empty string

*Evidence: Wikipedia encoding table, Cock et al. (2009)*

#### M3. Quality Score Decoding
- **M3.1** DecodeQualityScores Phred33: '!' (ASCII 33) → Q0
- **M3.2** DecodeQualityScores Phred33: 'I' (ASCII 73) → Q40
- **M3.3** DecodeQualityScores Phred64: '@' (ASCII 64) → Q0
- **M3.4** DecodeQualityScores Phred64: 'h' (ASCII 104) → Q40
- **M3.5** Empty quality string returns empty array

*Evidence: Wikipedia ASCII encoding table*

#### M4. Phred Mathematics
- **M4.1** PhredToErrorProbability(10) ≈ 0.1
- **M4.2** PhredToErrorProbability(20) ≈ 0.01
- **M4.3** PhredToErrorProbability(30) ≈ 0.001
- **M4.4** ErrorProbabilityToPhred(0.1) = 10
- **M4.5** ErrorProbabilityToPhred(0.01) = 20
- **M4.6** ErrorProbabilityToPhred(0.001) = 30

*Evidence: Wikipedia Phred score formula Q = -10×log₁₀(p)*

#### M5. Quality Score Encoding (Round-Trip)
- **M5.1** EncodeQualityScores Phred33: Q0 → '!', Q40 → 'I'
- **M5.2** EncodeQualityScores Phred64: Q0 → '@', Q40 → 'h'
- **M5.3** Decode(Encode(scores, enc), enc) == scores (round-trip)

*Evidence: Mathematical inverse property*

#### M6. Filtering
- **M6.1** FilterByQuality removes records below threshold
- **M6.2** FilterByQuality keeps records at/above threshold
- **M6.3** FilterByLength(min) removes short sequences
- **M6.4** FilterByLength(min, max) removes sequences outside range
- **M6.5** FilterByQuality on empty input returns empty

*Evidence: Implementation contract*

#### M7. Trimming
- **M7.1** TrimByQuality removes low-quality bases from 5' end
- **M7.2** TrimByQuality removes low-quality bases from 3' end
- **M7.3** TrimByQuality on all-high-quality returns unchanged
- **M7.4** TrimByQuality on all-low-quality returns empty sequence
- **M7.5** TrimAdapter removes adapter when found
- **M7.6** TrimAdapter returns unchanged when adapter not found

*Evidence: Standard quality trimming behavior*

#### M8. Statistics
- **M8.1** CalculateStatistics.TotalReads equals record count
- **M8.2** CalculateStatistics.TotalBases equals sum of sequence lengths
- **M8.3** CalculateStatistics.MeanReadLength = TotalBases / TotalReads
- **M8.4** CalculateStatistics.Q20Percentage in [0, 100]
- **M8.5** CalculateStatistics.Q30Percentage in [0, 100]
- **M8.6** CalculateStatistics.GcContent in [0, 1]
- **M8.7** Empty input returns zero statistics

*Evidence: Statistical definitions*

#### M9. Paired-End Support
- **M9.1** InterleavePairedReads alternates R1, R2, R1, R2...
- **M9.2** SplitInterleavedReads separates into equal R1/R2 lists
- **M9.3** InterleavePairedReads on unequal lengths stops at shorter

*Evidence: Wikipedia paired-end format description*

#### M10. File I/O
- **M10.1** ParseFile on non-existent file returns empty (no exception)
- **M10.2** ParseFile on valid file returns correct records
- **M10.3** WriteToFile creates valid FASTQ format
- **M10.4** Parse(WriteToFile content) reproduces records (round-trip)

*Evidence: Implementation contract*

---

### SHOULD Tests

#### S1. Header Parsing Details
- **S1.1** Header with space separates Id and Description
- **S1.2** Header without space: Description is empty
- **S1.3** Description can contain special characters

#### S2. Multi-line Handling
- **S2.1** Multi-line sequence assembled correctly
- **S2.2** Multi-line quality assembled correctly
- **S2.3** Blank lines between records skipped

#### S3. Position Quality Statistics
- **S3.1** CalculatePositionQuality returns per-position stats
- **S3.2** Position numbering is 1-based
- **S3.3** Standard deviation calculated correctly

---

### COULD Tests

#### C1. Edge Cases
- **C1.1** Sequence containing '+' character parsed correctly
- **C1.2** Very long sequences (10kb+) handled
- **C1.3** Unicode characters in header handled gracefully

#### C2. Performance
- **C2.1** Large file parsing completes in reasonable time
- **C2.2** Memory usage scales linearly with file size

---

## Audit Summary

### Existing Tests (FastqParserTests.cs)
| Test | Coverage | Status |
|------|----------|--------|
| Parse_SimpleFastq_ReturnsCorrectRecords | M1.1, M1.2 | ✓ Adequate |
| Parse_EmptyContent_ReturnsEmpty | M1.4 | ✓ Adequate |
| Parse_NullContent_ReturnsEmpty | M1.5 | ✓ Adequate |
| Parse_WithNoDescription_ParsesCorrectly | S1.2 | ✓ Adequate |
| Parse_RecordSequenceLength_MatchesQualityLength | M1.3 | ✓ Adequate |
| DetectEncoding_Phred33_ReturnsPhred33 | M2.1 | ✓ Adequate |
| DetectEncoding_Phred64_ReturnsPhred64 | M2.2 | ✓ Adequate |
| DecodeQualityScores_Phred33_ReturnsCorrectScores | M3.1, M3.2 | ✓ Adequate |
| DecodeQualityScores_Phred64_ReturnsCorrectScores | M3.3, M3.4 | ✓ Adequate |
| FilterByQuality_FiltersLowQuality | M6.1 | ✓ Adequate |
| FilterByLength_FiltersShortReads | M6.3 | ✓ Adequate |
| FilterByLength_WithMaxLength_FiltersBoth | M6.4 | ✓ Adequate |
| TrimByQuality_TrimsLowQualityEnds | M7.1, M7.2 | ✓ Adequate |
| TrimAdapter_RemovesAdapter | M7.5 | ✓ Adequate |
| TrimAdapter_NoAdapter_ReturnsUnchanged | M7.6 | ✓ Adequate |
| CalculateStatistics_ReturnsCorrectStats | M8.1-M8.3 | ✓ Adequate |
| CalculateStatistics_VariousLengths_CorrectMinMax | M8 | ✓ Adequate |
| CalculatePositionQuality_ReturnsQualityPerPosition | S3.1 | ✓ Adequate |
| InterleavePairedReads_CombinesReads | M9.1 | ✓ Adequate |
| SplitInterleavedReads_SeparatesReads | M9.2 | ✓ Adequate |
| Parse_MultiplePlusLines_ParsesCorrectly | C1.1 | ✓ Adequate |
| Parse_EmptyRecords_Skipped | S2.3 | ✓ Adequate |
| FilterByQuality_EmptyInput_ReturnsEmpty | M6.5 | ✓ Adequate |
| ParseFile_NonexistentFile_ReturnsEmpty | M10.1 | ✓ Adequate |
| ParseFile_ValidFile_ParsesRecords | M10.2 | ✓ Adequate |

### Missing Tests (All Closed)
| Test ID | Description | Status |
|---------|-------------|--------|
| M2.3 | DetectEncoding defaults Phred33 for ambiguous range | ✅ Covered |
| M2.4 | DetectEncoding empty string returns Phred33 | ✅ Covered |
| M3.5 | Empty quality string returns empty array | ✅ Covered |
| M4.1-M4.6 | Phred↔ErrorProbability conversion tests | ✅ Covered |
| M5.1-M5.3 | EncodeQualityScores and round-trip | ✅ Covered |
| M6.2 | FilterByQuality keeps records at threshold | ✅ Covered |
| M7.3-M7.4 | TrimByQuality edge cases | ✅ Covered |
| M8.4-M8.7 | Statistics invariants and edge cases | ✅ Covered |
| M9.3 | InterleavePairedReads unequal lengths | ✅ Covered |
| M10.3-M10.4 | Write and round-trip tests | ✅ Covered |
| S1.1, S1.3 | Header parsing edge cases | ✅ Covered |

### MCP Wrapper Tests (Seqeron.Mcp.Parsers.Tests)
Smoke tests exist for MCP bindings - no duplication with canonical tests.

---

## Consolidation Plan

1. **Canonical file:** `FastqParserTests.cs` - all deep logic tests
2. **MCP tests:** Keep as smoke tests only (schema validation + basic binding)
3. **No duplicates** between canonical and MCP test files
4. **Add missing tests** to canonical file

---

## Test Data

### Standard Test Record (Phred+33)
```csharp
const string StandardFastq = @"@SEQ_ID description
GATCGATCGATCGATC
+
IIIIIIIIIIIIIIII";
```

### Multi-Quality Test Records
```csharp
const string MultiQualityFastq = @"@high_quality
ACGT
+
IIII
@low_quality
ACGT
+
!!!!";
```

### Phred+64 Test Record
```csharp
const string Phred64Fastq = @"@read1
ACGT
+
hhhh";
```
