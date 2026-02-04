# PARSE-BED-001: BED File Parsing — Test Specification

**Test Unit ID:** PARSE-BED-001  
**Created:** 2026-02-05  
**Algorithm Group:** FileIO  
**Canonical Class:** BedParser

---

## 1. Scope

### 1.1 Canonical Methods (Deep Testing)
| Method | Type | Test Depth |
|--------|------|------------|
| `Parse(content, format)` | Canonical | Full |
| `ParseFile(filePath, format)` | File I/O | Full |
| `FilterByChrom(records, chrom)` | Filter | Full |
| `FilterByRegion(records, chrom, start, end)` | Filter | Full |
| `MergeOverlapping(records)` | Interval Op | Full |
| `Intersect(recordsA, recordsB)` | Interval Op | Full |

### 1.2 Supporting Methods (Standard Testing)
| Method | Type | Test Depth |
|--------|------|------------|
| `FilterByStrand(records, strand)` | Filter | Standard |
| `FilterByLength(records, min, max)` | Filter | Standard |
| `FilterByScore(records, min, max)` | Filter | Standard |
| `Subtract(recordsA, recordsB)` | Interval Op | Standard |
| `ExpandIntervals(records, upstream, downstream)` | Transform | Standard |
| `ExpandBlocks(record)` | BED12 | Standard |
| `GetIntrons(record)` | BED12 | Standard |
| `GetTotalBlockLength(record)` | BED12 | Standard |
| `CalculateStatistics(records)` | Stats | Standard |
| `WriteToStream(writer, records, format)` | Output | Standard |
| `Sort(records)` | Utility | Standard |

### 1.3 Wrapper/Delegate Methods (Smoke Testing)
| Method | Location | Test Depth |
|--------|----------|------------|
| `ParsersTools.BedParse` | MCP Parsers | Smoke (1-2 tests) |
| `ParsersTools.BedFilter` | MCP Parsers | Smoke (1-2 tests) |
| `ParsersTools.BedMerge` | MCP Parsers | Smoke (1-2 tests) |
| `ParsersTools.BedIntersect` | MCP Parsers | Smoke (1-2 tests) |
| `SequenceIO.ParseBedString` | SequenceIO | Smoke (existing) |

---

## 2. Must Tests (Evidence-Based)

### 2.1 Coordinate System Validation
| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M1 | Coordinate_ZeroBased_FirstBaseIsZero | UCSC FAQ | "The first base in a chromosome is numbered 0" |
| M2 | Coordinate_EndIsNonInclusive_LengthCorrect | UCSC FAQ | "chromEnd base is not included in the display" |
| M3 | Coordinate_ZeroLength_ValidForInsertions | UCSC FAQ | "chromStart=0, chromEnd=0 to represent an insertion" |
| M4 | Coordinate_First100Bases_StartZeroEndHundred | Wikipedia | "chr1:1-100 in browser = chromStart=0, chromEnd=100" |

### 2.2 Format Variants
| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M5 | Parse_BED3_ThreeColumnsRequired | UCSC FAQ | "three required fields" |
| M6 | Parse_BED6_IncludesNameScoreStrand | UCSC FAQ | Fields 4-6 definition |
| M7 | Parse_BED12_IncludesBlocks | UCSC FAQ | Fields 10-12 definition |
| M8 | Parse_MixedFormatsNotAllowed_ConsistentColumns | Wikipedia | "Each row must have same number of columns" |

### 2.3 Score Validation
| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M9 | Parse_Score_ClampedToRange0To1000 | UCSC FAQ | "A score between 0 and 1000" |
| M10 | Parse_Score_HigherThan1000_ClampedTo1000 | UCSC FAQ | Score must be in valid range |

### 2.4 Strand Validation
| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M11 | Parse_Strand_Plus_Valid | UCSC FAQ | "+" is valid strand |
| M12 | Parse_Strand_Minus_Valid | UCSC FAQ | "-" is valid strand |
| M13 | Parse_Strand_Dot_NoStrand | UCSC FAQ | "." indicates no strand |

### 2.5 Header Line Handling
| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M14 | Parse_TrackLine_Skipped | UCSC FAQ | "lines begin with the word 'browser' or 'track'" |
| M15 | Parse_BrowserLine_Skipped | UCSC FAQ | Browser settings lines |
| M16 | Parse_CommentLine_Skipped | Wikipedia | "#" denotes comment |

### 2.6 BED12 Block Validation
| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M17 | Parse_BED12_FirstBlockStartMustBeZero | UCSC FAQ | "first blockStart value must be 0" |
| M18 | Parse_BED12_BlocksNotOverlap | UCSC FAQ | "Blocks may not overlap" |
| M19 | Parse_BED12_BlockCountMatchesArrays | UCSC FAQ | "number of items... should correspond to blockCount" |
| M20 | Parse_BED12_FinalBlockMustReachEnd | UCSC FAQ | "final blockStart + final blockSize must equal chromEnd" |

### 2.7 Interval Operations
| # | Test | Evidence | Rationale |
|---|------|----------|-----------|
| M21 | MergeOverlapping_AdjacentIntervals_Merged | BEDTools | Standard merge behavior |
| M22 | MergeOverlapping_OverlappingIntervals_Merged | BEDTools | Overlapping regions combine |
| M23 | MergeOverlapping_DifferentChromosomes_NotMerged | BEDTools | Only same-chrom merges |
| M24 | Intersect_OverlappingRegions_ReturnsIntersection | BEDTools | Intersection semantics |
| M25 | Intersect_NoOverlap_ReturnsEmpty | BEDTools | No false positives |

---

## 3. Should Tests (Implementation Quality)

### 3.1 Delimiter Handling
| # | Test | Rationale |
|---|------|-----------|
| S1 | Parse_TabSeparated_ParsesCorrectly | Standard delimiter |
| S2 | Parse_SpaceSeparated_ParsesCorrectly | Allowed per UCSC |
| S3 | Parse_MixedDelimiters_HandleGracefully | Robustness |

### 3.2 Edge Cases
| # | Test | Rationale |
|---|------|-----------|
| S4 | Parse_EmptyContent_ReturnsEmpty | Empty input handling |
| S5 | Parse_NullContent_ReturnsEmpty | Null safety |
| S6 | Parse_InvalidCoordinates_SkipsLine | Malformed data handling |
| S7 | Parse_TooFewColumns_SkipsLine | Insufficient data |
| S8 | Parse_WhitespaceOnly_ReturnsEmpty | Empty lines |

### 3.3 Filter Operations
| # | Test | Rationale |
|---|------|-----------|
| S9 | FilterByChrom_CaseInsensitive | Chromosome name robustness |
| S10 | FilterByRegion_PartialOverlap_Included | Overlap semantics |
| S11 | FilterByStrand_InvalidStrand_NoMatch | Invalid input handling |

### 3.4 Statistics
| # | Test | Rationale |
|---|------|-----------|
| S12 | CalculateStatistics_EmptyRecords_ZeroCounts | Empty input |
| S13 | CalculateStatistics_ChromosomeCounts_Correct | Per-chrom counting |

### 3.5 I/O Operations
| # | Test | Rationale |
|---|------|-----------|
| S14 | ParseFile_NonexistentFile_ReturnsEmpty | File not found |
| S15 | ParseFile_ValidFile_ParsesRecords | File I/O |
| S16 | WriteAndRead_Roundtrip_PreservesData | Data integrity |

---

## 4. Could Tests (Enhancement)

### 4.1 Performance
| # | Test | Rationale |
|---|------|-----------|
| C1 | Parse_LargeFile_CompletesInReasonableTime | Performance baseline |

### 4.2 Extended Validation
| # | Test | Rationale |
|---|------|-----------|
| C2 | Parse_BED12_ExonCoordinatesCalculatedCorrectly | Block → absolute coords |
| C3 | GetIntrons_MultipleBlocks_IntronRegionsCorrect | Intron extraction |
| C4 | ExpandIntervals_StrandAware_UpstreamDownstream | Strand-aware expansion |

---

## 5. Test Audit (Existing Coverage)

### 5.1 BedParserTests.cs (Canonical - 42 tests)
| Category | Tests | Coverage | Status |
|----------|-------|----------|--------|
| Basic Parsing | 7 | BED3/6/12, empty, null, headers, comments | ✓ Covered |
| Record Properties | 4 | Length, HasBlocks, Score clamping | ✓ Covered |
| Filtering | 6 | Chrom, region, strand, length, score | ✓ Covered |
| Interval Operations | 6 | ToIntervals, Overlaps, Intersect, Merge, Subtract, Expand | ✓ Covered |
| Block Operations | 3 | ExpandBlocks, TotalBlockLength, GetIntrons | ✓ Covered |
| Statistics | 3 | Basic stats, chrom counts, empty | ✓ Covered |
| Writing | 2 | WriteToStream, Roundtrip | ✓ Covered |
| Edge Cases | 3 | Space-separated, invalid coords, few columns | ✓ Covered |
| File I/O | 2 | Nonexistent, valid file | ✓ Covered |
| Utility | 3 | Sort, Coverage, ExtractSequence | ✓ Covered |

### 5.2 MCP Parsers Tests (Delegate - Smoke)
| File | Tests | Status |
|------|-------|--------|
| BedParseTests.cs | 5 | ✓ Adequate smoke coverage |
| BedFilterTests.cs | 8 | ✓ Adequate smoke coverage |
| BedMergeTests.cs | 4 | ✓ Adequate smoke coverage |
| BedIntersectTests.cs | 5 | ✓ Adequate smoke coverage |

### 5.3 SequenceIO Tests (Legacy Wrapper)
| Tests | Status |
|-------|--------|
| 5 | ✓ Adequate for wrapper |

---

## 6. Consolidation Plan

### 6.1 Missing Tests to Add (Canonical)
| # | Test | Justification |
|---|------|---------------|
| 1 | Coordinate_ZeroBased_FirstBaseIsZero | M1 - explicit coordinate system test |
| 2 | Coordinate_EndIsNonInclusive_LengthCorrect | M2 - verify length = end - start |
| 3 | Coordinate_ZeroLength_ValidForInsertions | M3 - zero-length feature test |
| 4 | Parse_BED12_FirstBlockStartMustBeZero | M17 - block constraint |
| 5 | Parse_BED12_FinalBlockMustReachEnd | M20 - block constraint |
| 6 | Parse_BED12_BlockCountMatchesArrays | M19 - validate arrays |
| 7 | GenomicInterval_Union_ReturnsUnion | Complete interval API |
| 8 | ExpandIntervals_ZeroExpansion_Unchanged | Edge case |

### 6.2 Tests to Refactor
- None needed - existing structure is clean

### 6.3 Duplicates to Remove
- None identified - canonical vs MCP tests have clear separation

---

## 7. Decisions

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Score values >1000 are clamped to 1000 | Per UCSC spec, implementation already does this |
| D2 | Invalid coordinate lines are skipped | Robustness over strict parsing |
| D3 | Both tab and space delimiters accepted | Per UCSC allowing both |
| D4 | Chromosome comparison is case-insensitive | Common practice in genomics |

---

## 8. Canonical Test File

**Location:** `tests/Seqeron/Seqeron.Genomics.Tests/BedParserTests.cs`

All deep, evidence-based tests for BedParser consolidated here.
MCP wrapper tests remain separate with smoke-level coverage only.
