# bed_intersect

Find intersecting regions between two BED datasets.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `bed_intersect` |
| **Method ID** | `BedParser.Intersect` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds intersecting (overlapping) regions between two sets of BED records. Returns the overlapping portions of features from set A that overlap with features in set B. Similar to bedtools intersect.

## Core Documentation Reference

- Source: [BedParser.cs#L319](../../../../Seqeron.Genomics/BedParser.cs#L319)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `contentA` | string | Yes | First BED format content (features to intersect) |
| `contentB` | string | Yes | Second BED format content (reference features) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of intersection BED records |
| `intersectionCount` | integer | Number of intersections found |
| `countA` | integer | Number of records in set A |
| `countB` | integer | Number of records in set B |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content A cannot be null or empty |
| 1002 | Content B cannot be null or empty |

## Examples

### Example 1: Find overlapping regions

**User Prompt:**
> Find where these two sets of features overlap

**Expected Tool Call:**
```json
{
  "tool": "bed_intersect",
  "arguments": {
    "contentA": "chr1\t100\t200",
    "contentB": "chr1\t150\t250"
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 150, "chromEnd": 200, "length": 50 }
  ],
  "intersectionCount": 1,
  "countA": 1,
  "countB": 1
}
```

### Example 2: No overlap

**Expected Tool Call:**
```json
{
  "tool": "bed_intersect",
  "arguments": {
    "contentA": "chr1\t100\t200",
    "contentB": "chr1\t300\t400"
  }
}
```

**Response:**
```json
{
  "records": [],
  "intersectionCount": 0,
  "countA": 1,
  "countB": 1
}
```

### Example 3: Multiple overlaps

**Expected Tool Call:**
```json
{
  "tool": "bed_intersect",
  "arguments": {
    "contentA": "chr1\t100\t300",
    "contentB": "chr1\t150\t200\nchr1\t250\t350"
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 150, "chromEnd": 200, "length": 50 },
    { "chrom": "chr1", "chromStart": 250, "chromEnd": 300, "length": 50 }
  ],
  "intersectionCount": 2,
  "countA": 1,
  "countB": 2
}
```

## Performance

- **Time Complexity:** O(n × m) in worst case, but optimized with chromosome indexing
- **Space Complexity:** O(n + m) for storing records and results

## See Also

- [bed_parse](bed_parse.md) - Parse BED format
- [bed_filter](bed_filter.md) - Filter BED records
- [bed_merge](bed_merge.md) - Merge overlapping records
