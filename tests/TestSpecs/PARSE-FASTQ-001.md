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
- **M3.3** DecodeQualityScores Phred33: '~' (ASCII 126) → Q93
- **M3.4** DecodeQualityScores Phred64: '@' (ASCII 64) → Q0
- **M3.5** DecodeQualityScores Phred64: 'h' (ASCII 104) → Q40
- **M3.6** DecodeQualityScores Phred64: '~' (ASCII 126) → Q62
- **M3.7** Empty quality string returns empty array

*Evidence: Wikipedia FASTQ encoding table — Sanger ASCII 33-126 (Q 0-93), Phred+64 ASCII 64-126 (Q 0-62)*

#### M4. Phred Mathematics
- **M4.1** PhredToErrorProbability(0) ≈ 1.0
- **M4.2** PhredToErrorProbability(10) ≈ 0.1
- **M4.3** PhredToErrorProbability(20) ≈ 0.01
- **M4.4** PhredToErrorProbability(30) ≈ 0.001
- **M4.5** PhredToErrorProbability(40) ≈ 0.0001
- **M4.6** ErrorProbabilityToPhred(0.1) = 10
- **M4.7** ErrorProbabilityToPhred(0.01) = 20
- **M4.8** ErrorProbabilityToPhred(0.001) = 30
- **M4.9** ErrorProbabilityToPhred(0) = 93 (max Sanger representable)

*Evidence: Wikipedia Phred quality score Q = -10×log₁₀(p); Symbols table Q0=1.000…Q40=0.0001; max Q93 per Sanger ASCII 126-33*

#### M5. Quality Score Encoding (Round-Trip)
- **M5.1** EncodeQualityScores Phred33: Q0 → '!', Q40 → 'I', Q93 → '~'
- **M5.2** EncodeQualityScores Phred64: Q0 → '@', Q40 → 'h', Q62 → '~'
- **M5.3** Decode(Encode(scores, Phred33), Phred33) == scores for range [0, 93]
- **M5.4** Decode(Encode(scores, Phred64), Phred64) == scores for range [0, 62]

*Evidence: Wikipedia FASTQ encoding — Sanger Q 0-93 (ASCII 33-126), Phred+64 Q 0-62 (ASCII 64-126)*

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

## Deviations and Assumptions

None. All behavior is evidence-backed:

| # | Behavior | Justification | Source |
|---|----------|---------------|--------|
| 1 | Detection heuristic: chars < '@' → Phred+33; chars > 'I' → Phred+64; default Phred+33 | Standard auto-detection approach | Wikipedia FASTQ Encoding |
| 2 | Phred+33 encodes Q 0-93 (ASCII 33-126) | Sanger format full range; PacBio HiFi uses up to Q93 | Wikipedia FASTQ Encoding chart |
| 3 | Phred+64 encodes Q 0-62 (ASCII 64-126) | Illumina 1.3-1.7 full range | Wikipedia FASTQ Encoding chart |
| 4 | ErrorProbabilityToPhred(0) = 93 | Q = -10×log₁₀(0) = ∞; capped at max representable Sanger value | Wikipedia FASTQ: Sanger ASCII 33-126 |
| 5 | '+' allowed within sequence data | Parser reads sequence until standalone '+' line | Wikipedia FASTQ Format |
| 6 | Multi-line sequence/quality supported | Reads until '+' / sequence-length reached | Wikipedia FASTQ Format |

---

## Audit Summary

### Tests (FastqParserTests.cs — 60 tests)
| Test | Coverage |
|------|----------|
| Parse_SimpleFastq_ReturnsCorrectRecords | M1.1, M1.2 |
| Parse_EmptyContent_ReturnsEmpty | M1.4 |
| Parse_NullContent_ReturnsEmpty | M1.5 |
| Parse_WithNoDescription_ParsesCorrectly | S1.2 |
| Parse_RecordSequenceLength_MatchesQualityLength | M1.3 |
| Parse_HeaderWithSpace_SeparatesIdAndDescription | S1.1 |
| Parse_DescriptionWithSpecialCharacters_ParsedCorrectly | S1.3 |
| Parse_MultiLineSequence_AssembledCorrectly | S2.1 |
| Parse_MultiLineQuality_AssembledCorrectly | S2.2 |
| DetectEncoding_Phred33_ReturnsPhred33 | M2.1 |
| DetectEncoding_Phred64_ReturnsPhred64 | M2.2 |
| DetectEncoding_AmbiguousRange_DefaultsToPhred33 | M2.3 |
| DetectEncoding_EmptyString_ReturnsPhred33 | M2.4 |
| DecodeQualityScores_Phred33_ReturnsCorrectScores | M3.1-M3.3 |
| DecodeQualityScores_Phred64_ReturnsCorrectScores | M3.4-M3.6 |
| DecodeQualityScores_EmptyString_ReturnsEmptyArray | M3.7 |
| DecodeQualityScores_NullString_ReturnsEmptyArray | M3.7 |
| PhredToErrorProbability_Q0_Returns1 | M4.1 |
| PhredToErrorProbability_Q10_Returns0Point1 | M4.2 |
| PhredToErrorProbability_Q20_Returns0Point01 | M4.3 |
| PhredToErrorProbability_Q30_Returns0Point001 | M4.4 |
| PhredToErrorProbability_Q40_Returns0Point0001 | M4.5 |
| ErrorProbabilityToPhred_0Point1_ReturnsQ10 | M4.6 |
| ErrorProbabilityToPhred_0Point01_ReturnsQ20 | M4.7 |
| ErrorProbabilityToPhred_0Point001_ReturnsQ30 | M4.8 |
| ErrorProbabilityToPhred_ZeroOrNegative_ReturnsMaxQuality | M4.9 |
| EncodeQualityScores_Phred33_EncodesCorrectly | M5.1 |
| EncodeQualityScores_Phred64_EncodesCorrectly | M5.2 |
| EncodeDecodeRoundTrip_Phred33_PreservesScores | M5.3 |
| EncodeDecodeRoundTrip_Phred64_PreservesScores | M5.4 |
| FilterByQuality_FiltersLowQuality | M6.1, M6.2 |
| FilterByQuality_KeepsRecordsAtThreshold | M6.2 |
| FilterByQuality_EmptyInput_ReturnsEmpty | M6.5 |
| FilterByLength_FiltersShortReads | M6.3 |
| FilterByLength_WithMaxLength_FiltersBoth | M6.4 |
| TrimByQuality_TrimsLowQualityEnds | M7.1, M7.2 |
| TrimByQuality_AllHighQuality_ReturnsUnchanged | M7.3 |
| TrimByQuality_AllLowQuality_ReturnsEmptySequence | M7.4 |
| TrimAdapter_RemovesAdapter | M7.5 |
| TrimAdapter_NoAdapter_ReturnsUnchanged | M7.6 |
| CalculateStatistics_ReturnsCorrectStats | M8.1-M8.3 |
| CalculateStatistics_VariousLengths_CorrectMinMax | M8 |
| CalculateStatistics_EmptyInput_ReturnsZeros | M8.7 |
| CalculateStatistics_Q20Percentage_InValidRange | M8.4 |
| CalculateStatistics_Q30Percentage_InValidRange | M8.5 |
| CalculateStatistics_GcContent_InValidRange | M8.6 |
| CalculateStatistics_HighQualityReads_HasHighQ30 | M8.5 |
| CalculatePositionQuality_ReturnsQualityPerPosition | S3.1, S3.2, S3.3 |
| InterleavePairedReads_CombinesReads | M9.1 |
| SplitInterleavedReads_SeparatesReads | M9.2 |
| InterleavePairedReads_UnequalLengths_StopsAtShorter | M9.3 |
| Parse_MultiplePlusLines_ParsesCorrectly | C1.1 |
| Parse_VeryLongSequence_HandledCorrectly | C1.2 |
| Parse_UnicodeInHeader_HandledGracefully | C1.3 |
| Parse_EmptyRecords_Skipped | S2.3 |
| ParseFile_NonexistentFile_ReturnsEmpty | M10.1 |
| ParseFile_ValidFile_ParsesRecords | M10.2 |
| WriteToFile_CreatesValidFastq | M10.3 |
| WriteAndParseRoundTrip_PreservesRecords | M10.4 |
| ToFastqString_FormatsCorrectly | M10 |

### MCP Wrapper Tests (Seqeron.Mcp.Parsers.Tests — 36 tests)
Smoke tests for MCP bindings — no duplication with canonical tests.

---

### Actions Taken

| # | Action | Test | Detail |
|---|--------|------|--------|
| 1 | ⚠ Strengthened | FilterByQuality_FiltersLowQuality | exact count=2 + exact IDs instead of `Count < records.Count` |
| 2 | ⚠ Strengthened | TrimByQuality_TrimsLowQualityEnds | exact sequence "GTACGTAC" + exact quality instead of `Length < 12` |
| 3 | ⚠ Strengthened | TrimAdapter_RemovesAdapter | exact remaining "ACGTACGTACGTAAA" instead of `Does.Not.Contain` |
| 4 | ⚠ Strengthened | CalculateStatistics_Q20Percentage_InValidRange | exact 100.0 instead of `InRange(0, 100)` |
| 5 | ⚠ Strengthened | CalculateStatistics_Q30Percentage_InValidRange | exact 100.0 instead of `InRange(0, 100)` |
| 6 | ⚠ Strengthened | CalculateStatistics_GcContent_InValidRange | exact 0.5 instead of `InRange(0, 1)` |
| 7 | ⚠ Strengthened | SplitInterleavedReads_SeparatesReads | verify exact sequences per read pair |
| 8 | ⚠ Strengthened | ParseFile_ValidFile_ParsesRecords | verify record Id+Sequence, not just count |
| 9 | ⚠ Strengthened | WriteToFile_CreatesValidFastq | exact line-by-line ReadAllLines check |
| 10 | ⚠ Strengthened | ToFastqString_FormatsCorrectly | exact 4-line split verification |
| 11 | ⚠ Strengthened | CalculatePositionQuality_ReturnsQualityPerPosition | exact mean=39.5, stddev=0.5, 1-based positions |
| 12 | ❌ Added | Parse_DescriptionWithSpecialCharacters_ParsedCorrectly | S1.3 |
| 13 | ❌ Added | Parse_MultiLineSequence_AssembledCorrectly | S2.1 |
| 14 | ❌ Added | Parse_MultiLineQuality_AssembledCorrectly | S2.2 |
| 15 | ❌ Added | Parse_VeryLongSequence_HandledCorrectly | C1.2 |
| 16 | ❌ Added | Parse_UnicodeInHeader_HandledGracefully | C1.3 |
| 17 | 🔁 Removed | FilterByQuality_RemovesRecordsBelowThreshold | duplicate of M6.1 (covered by strengthened FiltersLowQuality) |
