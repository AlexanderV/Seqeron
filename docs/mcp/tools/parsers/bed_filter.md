# bed_filter

Filter BED records by chromosome, region, strand, length, or score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `bed_filter` |
| **Method ID** | `BedParser.FilterBy*` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Filters BED records using multiple optional criteria. All filters can be combined to create complex queries. Returns records matching all specified criteria.

**Available Filters:**
- **Chromosome**: Filter by chromosome name
- **Region**: Filter by genomic coordinates (overlapping)
- **Strand**: Filter by + or - strand
- **Length**: Filter by minimum/maximum feature length
- **Score**: Filter by minimum/maximum score

## Core Documentation Reference

- Source: [BedParser.cs#L217](../../../../src/Seqeron/Seqeron.Genomics/BedParser.cs#L217)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | BED format content to filter |
| `chrom` | string | No | - | Filter by chromosome name |
| `regionStart` | integer | No | - | Region start (requires chrom and regionEnd) |
| `regionEnd` | integer | No | - | Region end (requires chrom and regionStart) |
| `strand` | string | No | - | Filter by strand: "+" or "-" |
| `minLength` | integer | No | - | Minimum feature length |
| `maxLength` | integer | No | - | Maximum feature length |
| `minScore` | integer | No | - | Minimum score |
| `maxScore` | integer | No | - | Maximum score |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `records` | array | List of filtered BED records |
| `passedCount` | integer | Number of records that passed filters |
| `totalCount` | integer | Total number of input records |
| `passedPercentage` | number | Percentage of records that passed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Strand must be '+' or '-' |

## Examples

### Example 1: Filter by chromosome

**User Prompt:**
> Show me all features on chromosome 1

**Expected Tool Call:**
```json
{
  "tool": "bed_filter",
  "arguments": {
    "content": "chr1\t100\t200\nchr2\t300\t400\nchr1\t500\t600",
    "chrom": "chr1"
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 100, "chromEnd": 200, "length": 100 },
    { "chrom": "chr1", "chromStart": 500, "chromEnd": 600, "length": 100 }
  ],
  "passedCount": 2,
  "totalCount": 3,
  "passedPercentage": 66.67
}
```

### Example 2: Filter by region

**User Prompt:**
> Find features overlapping chr1:150-350

**Expected Tool Call:**
```json
{
  "tool": "bed_filter",
  "arguments": {
    "content": "chr1\t100\t200\tA\nchr1\t300\t400\tB\nchr1\t500\t600\tC",
    "chrom": "chr1",
    "regionStart": 150,
    "regionEnd": 350
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 100, "chromEnd": 200, "length": 100, "name": "A" },
    { "chrom": "chr1", "chromStart": 300, "chromEnd": 400, "length": 100, "name": "B" }
  ],
  "passedCount": 2,
  "totalCount": 3,
  "passedPercentage": 66.67
}
```

### Example 3: Combine multiple filters

**User Prompt:**
> Find plus-strand features on chr1 with score >= 500

**Expected Tool Call:**
```json
{
  "tool": "bed_filter",
  "arguments": {
    "content": "chr1\t100\t200\tA\t600\t+\nchr1\t300\t400\tB\t400\t+\nchr1\t500\t600\tC\t700\t-",
    "chrom": "chr1",
    "strand": "+",
    "minScore": 500
  }
}
```

**Response:**
```json
{
  "records": [
    { "chrom": "chr1", "chromStart": 100, "chromEnd": 200, "length": 100, "name": "A", "score": 600, "strand": "+" }
  ],
  "passedCount": 1,
  "totalCount": 3,
  "passedPercentage": 33.33
}
```

## Performance

- **Time Complexity:** O(n) where n is number of records
- **Space Complexity:** O(n) for storing filtered records

## See Also

- [bed_parse](bed_parse.md) - Parse BED format
- [bed_merge](bed_merge.md) - Merge overlapping BED records
- [bed_intersect](bed_intersect.md) - Find intersecting regions
