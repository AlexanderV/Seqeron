# bed_merge

Merge overlapping BED records into single intervals.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `bed_merge` |
| **Method ID** | `BedParser.MergeOverlapping` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Merges overlapping or adjacent BED records on the same chromosome into single intervals. Records are sorted by chromosome and position before merging. Useful for consolidating feature annotations or reducing redundancy.

## Core Documentation Reference

- Source: [BedParser.cs#L284](../../../../Seqeron.Genomics/BedParser.cs#L284)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `content` | string | Yes | BED format content to merge |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of merged BED records |
| `mergedCount` | integer | Number of records after merging |
| `originalCount` | integer | Number of records before merging |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |

## Examples

### Example 1: Merge overlapping intervals

**User Prompt:**
> Merge these overlapping genomic regions

**Expected Tool Call:**
```json
{
  "tool": "bed_merge",
  "arguments": {
    "content": "chr1\t100\t200\nchr1\t150\t250\nchr1\t400\t500"
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 100, "chromEnd": 250, "length": 150 },
    { "chrom": "chr1", "chromStart": 400, "chromEnd": 500, "length": 100 }
  ],
  "mergedCount": 2,
  "originalCount": 3
}
```

### Example 2: Different chromosomes (no merge)

**Expected Tool Call:**
```json
{
  "tool": "bed_merge",
  "arguments": {
    "content": "chr1\t100\t200\nchr2\t100\t200"
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 100, "chromEnd": 200, "length": 100 },
    { "chrom": "chr2", "chromStart": 100, "chromEnd": 200, "length": 100 }
  ],
  "mergedCount": 2,
  "originalCount": 2
}
```

## Performance

- **Time Complexity:** O(n log n) for sorting + O(n) for merging
- **Space Complexity:** O(n) for storing sorted and merged records

## See Also

- [bed_parse](bed_parse.md) - Parse BED format
- [bed_filter](bed_filter.md) - Filter BED records
- [bed_intersect](bed_intersect.md) - Find intersecting regions
