# bed_parse

Parse BED format content into genomic region records.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `bed_parse` |
| **Method ID** | `BedParser.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses BED (Browser Extensible Data) format content into genomic region records. BED format is widely used to describe genomic features like genes, regulatory elements, and alignment regions. Supports BED3, BED6, and BED12 formats.

**BED Format Variants:**
- **BED3**: chrom, chromStart, chromEnd (minimal)
- **BED6**: + name, score, strand
- **BED12**: + thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts

## Core Documentation Reference

- Source: [BedParser.cs#L102](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/BedParser.cs#L102)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | BED format content to parse |
| `format` | string | No | `"auto"` | BED format: `"bed3"`, `"bed6"`, `"bed12"`, or `"auto"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of parsed BED records |
| `records[].chrom` | string | Chromosome name |
| `records[].chromStart` | integer | Start position (0-based) |
| `records[].chromEnd` | integer | End position (exclusive) |
| `records[].length` | integer | Feature length (chromEnd - chromStart) |
| `records[].name` | string? | Feature name (BED4+) |
| `records[].score` | integer? | Score 0-1000 (BED5+) |
| `records[].strand` | string? | Strand "+" or "-" (BED6+) |
| `records[].thickStart` | integer? | Thick drawing start (BED12) |
| `records[].thickEnd` | integer? | Thick drawing end (BED12) |
| `records[].itemRgb` | string? | RGB color (BED12) |
| `records[].blockCount` | integer? | Number of blocks (BED12) |
| `records[].blockSizes` | integer[]? | Block sizes (BED12) |
| `records[].blockStarts` | integer[]? | Block start positions (BED12) |
| `count` | integer | Number of records parsed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Invalid format: {format}. Use 'bed3', 'bed6', 'bed12', or 'auto' |

## Examples

### Example 1: Parse BED3 format

**User Prompt:**
> Parse this BED file with genomic regions

**Expected Tool Call:**
```json
{
  "tool": "bed_parse",
  "arguments": {
    "content": "chr1\t100\t200\nchr1\t300\t500"
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 100, "chromEnd": 200, "length": 100 },
    { "chrom": "chr1", "chromStart": 300, "chromEnd": 500, "length": 200 }
  ],
  "count": 2
}
```

### Example 2: Parse BED6 with gene names

**User Prompt:**
> Parse this BED6 file with gene annotations

**Expected Tool Call:**
```json
{
  "tool": "bed_parse",
  "arguments": {
    "content": "chr1\t1000\t5000\tGENE1\t900\t+\nchr1\t6000\t8000\tGENE2\t800\t-",
    "format": "bed6"
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 1000, "chromEnd": 5000, "length": 4000, "name": "GENE1", "score": 900, "strand": "+" },
    { "chrom": "chr1", "chromStart": 6000, "chromEnd": 8000, "length": 2000, "name": "GENE2", "score": 800, "strand": "-" }
  ],
  "count": 2
}
```

## Performance

- **Time Complexity:** O(n) where n is total content length
- **Space Complexity:** O(n) for storing parsed records

## Notes

- Comment lines (starting with #) and track lines are automatically skipped
- Coordinates are 0-based, half-open (start inclusive, end exclusive)
- Empty lines are ignored

## See Also

- [bed_filter](bed_filter.md) - Filter BED records by chromosome or region
- [bed_merge](bed_merge.md) - Merge overlapping BED records
- [bed_intersect](bed_intersect.md) - Find intersecting regions
