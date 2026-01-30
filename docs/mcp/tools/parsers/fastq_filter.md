# fastq_filter

Filter FASTQ reads by minimum average quality score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_filter` |
| **Method ID** | `FastqParser.FilterByQuality` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Filters FASTQ reads by minimum average quality score. Returns only reads that meet the quality threshold. Useful for quality control and removing low-quality reads before downstream analysis.

## Core Documentation Reference

- Source: [FastqParser.cs#L228](../../../../Seqeron.Genomics/FastqParser.cs#L228)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | FASTQ format content to filter |
| `minQuality` | number | Yes | - | Minimum average quality score threshold |
| `encoding` | string | No | `"auto"` | Quality encoding: `"phred33"`, `"phred64"`, or `"auto"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries` | array | List of filtered FASTQ records (same format as fastq_parse) |
| `passedCount` | integer | Number of reads that passed the filter |
| `totalCount` | integer | Total number of input reads |
| `passedPercentage` | number | Percentage of reads that passed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Minimum quality must be non-negative |
| 1003 | Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto' |

## Examples

### Example 1: Filter by Q30 threshold

**User Prompt:**
> Filter these reads to keep only those with average quality >= 30

**Expected Tool Call:**
```json
{
  "tool": "fastq_filter",
  "arguments": {
    "content": "@high_qual\nATGC\n+\nIIII\n@low_qual\nATGC\n+\n!!!!",
    "minQuality": 30,
    "encoding": "phred33"
  }
}
```

**Response:**
```json
{
  "entries": [
    {
      "id": "high_qual",
      "description": null,
      "sequence": "ATGC",
      "qualityString": "IIII",
      "qualityScores": [40, 40, 40, 40],
      "length": 4
    }
  ],
  "passedCount": 1,
  "totalCount": 2,
  "passedPercentage": 50.0
}
```

### Example 2: Filter with Q20 threshold (less stringent)

**User Prompt:**
> Keep reads with average quality at least 20

**Expected Tool Call:**
```json
{
  "tool": "fastq_filter",
  "arguments": {
    "content": "@read1\nACGT\n+\n5555\n@read2\nACGT\n+\n!!!!",
    "minQuality": 20
  }
}
```

**Response:**
```json
{
  "entries": [
    {
      "id": "read1",
      "description": null,
      "sequence": "ACGT",
      "qualityString": "5555",
      "qualityScores": [20, 20, 20, 20],
      "length": 4
    }
  ],
  "passedCount": 1,
  "totalCount": 2,
  "passedPercentage": 50.0
}
```

## Performance

- **Time Complexity:** O(n) where n is total content length
- **Space Complexity:** O(n) for storing parsed and filtered records

## Quality Thresholds

Common quality thresholds:
- **Q20**: 99% accuracy (1% error rate) - minimum for most analyses
- **Q30**: 99.9% accuracy (0.1% error rate) - high-quality threshold

## See Also

- [fastq_parse](fastq_parse.md) - Parse FASTQ into records
- [fastq_statistics](fastq_statistics.md) - Calculate quality statistics
