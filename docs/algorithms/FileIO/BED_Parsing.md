# BED Format Parsing

**Algorithm Group:** FileIO  
**Test Unit:** PARSE-BED-001

---

## Overview

BED (Browser Extensible Data) is a text-based format for storing genomic regions as coordinates with optional annotations. Developed during the Human Genome Project, it has become a de facto standard in bioinformatics.

---

## Format Specification

### Column Structure

| Column | Field | Required | Description |
|--------|-------|----------|-------------|
| 1 | chrom | Yes | Chromosome or scaffold name |
| 2 | chromStart | Yes | Start position (0-based) |
| 3 | chromEnd | Yes | End position (non-inclusive) |
| 4 | name | No | Feature name |
| 5 | score | No | Score 0-1000 |
| 6 | strand | No | +, -, or . |
| 7 | thickStart | No | CDS start |
| 8 | thickEnd | No | CDS end |
| 9 | itemRgb | No | RGB color |
| 10 | blockCount | No | Number of exons |
| 11 | blockSizes | No | Comma-separated block sizes |
| 12 | blockStarts | No | Comma-separated block starts |

### Coordinate System

BED uses a 0-based, half-open coordinate system:
- `chromStart`: 0-based (first base is position 0)
- `chromEnd`: Non-inclusive (the end base is NOT part of the feature)

**Example**: The first 100 bases of a chromosome:
- chromStart = 0, chromEnd = 100
- Represents bases at positions 0, 1, 2, ..., 99
- Browser display: chr1:1-100 (1-based for user display)

**Feature Length**: `chromEnd - chromStart`

### Zero-Length Features

When `chromStart == chromEnd`, the feature has zero length. This represents insertion points:
- `chr1 5 5` = insertion point at position 5

---

## BED Variants

| Format | Columns | Use Case |
|--------|---------|----------|
| BED3 | 3 | Simple intervals |
| BED4 | 4 | Named intervals |
| BED5 | 5 | Scored intervals |
| BED6 | 6 | Stranded intervals |
| BED12 | 12 | Gene models with exons |

---

## Header Lines

BED files may contain header lines (custom tracks only):

```
track name=myTrack description="My Track"
browser position chr7:127471196-127495720
# This is a comment
chr7	127471196	127472363	Feature1
```

Header line types:
- `track `: Display settings for genome browsers
- `browser `: Browser navigation settings
- `#`: Comments

**Note**: Header lines are NOT allowed in files processed by bedToBigBed.

---

## BED12 Block Structure

For gene models with exons:

```
chr1  1000  5000  gene1  900  +  1100  4900  0  3  100,200,300  0,1000,3700
```

**Block Constraints**:
1. First `blockStart` must be 0
2. `blockStarts` are relative to `chromStart`
3. Last `blockStart + blockSize` must equal `chromEnd - chromStart`
4. Blocks must not overlap

**Calculating Exon Coordinates**:
```
exon[i].start = chromStart + blockStarts[i]
exon[i].end = chromStart + blockStarts[i] + blockSizes[i]
```

---

## Implementation

### BedParser Class

```csharp
public static class BedParser
{
    // Parsing
    public static IEnumerable<BedRecord> Parse(string content, BedFormat format = BedFormat.Auto);
    public static IEnumerable<BedRecord> ParseFile(string filePath, BedFormat format = BedFormat.Auto);
    
    // Filtering
    public static IEnumerable<BedRecord> FilterByChrom(IEnumerable<BedRecord> records, string chrom);
    public static IEnumerable<BedRecord> FilterByRegion(IEnumerable<BedRecord> records, string chrom, int start, int end);
    public static IEnumerable<BedRecord> FilterByStrand(IEnumerable<BedRecord> records, char strand);
    
    // Interval Operations
    public static IEnumerable<BedRecord> MergeOverlapping(IEnumerable<BedRecord> records);
    public static IEnumerable<BedRecord> Intersect(IEnumerable<BedRecord> a, IEnumerable<BedRecord> b);
    public static IEnumerable<BedRecord> Subtract(IEnumerable<BedRecord> a, IEnumerable<BedRecord> b);
    
    // Block Operations (BED12)
    public static IEnumerable<BedRecord> ExpandBlocks(BedRecord record);
    public static int GetTotalBlockLength(BedRecord record);
    public static IEnumerable<BedRecord> GetIntrons(BedRecord record);
}
```

### BedRecord Structure

```csharp
public readonly record struct BedRecord(
    string Chrom,
    int ChromStart,
    int ChromEnd,
    string? Name = null,
    int? Score = null,
    char? Strand = null,
    int? ThickStart = null,
    int? ThickEnd = null,
    string? ItemRgb = null,
    int? BlockCount = null,
    int[]? BlockSizes = null,
    int[]? BlockStarts = null)
{
    public int Length => ChromEnd - ChromStart;
    public bool HasBlocks => BlockCount.HasValue && BlockSizes != null && BlockStarts != null;
}
```

---

## Complexity

| Operation | Time Complexity |
|-----------|-----------------|
| Parse | O(n) |
| FilterByChrom | O(n) |
| FilterByRegion | O(n) |
| MergeOverlapping | O(n log n) |
| Intersect | O(n × m) |
| Subtract | O(n × m) |

Where n = number of records, m = number of records in second set.

---

## References

1. UCSC Genome Browser FAQ - BED Format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
2. Wikipedia - BED (file format): https://en.wikipedia.org/wiki/BED_(file_format)
3. GA4GH BED v1.0 Specification (2021)
4. Kent WJ, et al. (2002). "The Human Genome Browser at UCSC". Genome Research.
