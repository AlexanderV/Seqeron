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

## 5. Test Audit (Current Coverage)

### 5.1 BedParserTests.cs (Canonical — 74 tests)
| Category | Tests | Coverage | Status |
|----------|-------|----------|--------|
| Basic Parsing | 8 | BED3/6/12, empty, null, whitespace-only, headers, comments | ✓ |
| Record Properties | 4 | Length, HasBlocks(true/false), Score clamping | ✓ |
| Filtering | 7 | Chrom, case-insensitive, region, strand, invalid strand, length, score | ✓ |
| Interval Operations | 10 | ToIntervals, Overlaps, Intersect (overlap + no-overlap), Merge (overlap + adjacent + diff chrom), Subtract, Expand | ✓ |
| Block Operations | 2 | TotalBlockLength, GetIntrons (with exact coordinates) | ✓ |
| Statistics | 3 | Basic stats, chrom counts, empty | ✓ |
| Writing | 2 | WriteToStream (exact line output), Roundtrip (all fields) | ✓ |
| Edge Cases | 6 | Space-delimited (with values), invalid coords, few columns, chromStart>End, mixed columns, negative score | ✓ |
| File I/O | 2 | Nonexistent, valid file | ✓ |
| Utility | 4 | Sort, Coverage (exact depths), ExtractSequence, ExtractSequence minus strand (non-palindrome RC) | ✓ |
| Coordinate System | 4 | Zero-based, non-inclusive end, zero-length, browser conversion | ✓ |
| BED12 Validation | 10 | First blockStart, blockCount match, final block, exon coords, UCSC examples, overlapping blocks | ✓ |
| GenomicInterval | 3 | Length, Union, non-overlap intersect | ✓ |
| Expand Intervals | 2 | Zero expansion, upstream clamp | ✓ |
| UCSC Reference Data | 3 | Gene structure, ChIP-seq style, hg38 coords | ✓ |
| Mutation-Killing | 2 | Subtract || → && mutations | ✓ |

### 5.2 MCP Parsers Tests (Delegate — Smoke)
| File | Tests | Status |
|------|-------|--------|
| BedParseTests.cs | 5 | ✓ |
| BedFilterTests.cs | 8 | ✓ |
| BedMergeTests.cs | 4 | ✓ |
| BedIntersectTests.cs | 5 | ✓ |

### 5.3 SequenceIO Tests (Legacy Wrapper)
| Tests | Status |
|-------|--------|
| 5 | ✓ |

---

## 6. Deviations and Assumptions

**None.** Implementation strictly follows UCSC Genome Browser BED specification and Wikipedia description.

All validation rules from the specification are enforced:
- **Coordinate validation**: chromStart ≤ chromEnd (UCSC FAQ)
- **Column consistency**: all data lines must have the same number of fields (UCSC FAQ + Wikipedia)
- **Score range**: clamped to [0, 1000] (UCSC FAQ)
- **Strand values**: only "+", "-", "." accepted (UCSC FAQ)
- **BED12 block constraints** (UCSC FAQ):
  - First blockStart must be 0
  - Final blockStart + final blockSize must equal chromEnd − chromStart
  - Blocks may not overlap
  - blockCount must match blockSizes/blockStarts array lengths
- **Header lines**: "track", "browser", "#" lines skipped (UCSC FAQ + Wikipedia)
- **Delimiters**: tab-separated preferred; space-separated accepted as fallback (UCSC FAQ + Wikipedia)

---

## 7. Canonical Test File

**Location:** `tests/Seqeron/Seqeron.Genomics.Tests/BedParserTests.cs`

All deep, evidence-based tests for BedParser consolidated here.
MCP wrapper tests remain separate with smoke-level coverage only.
