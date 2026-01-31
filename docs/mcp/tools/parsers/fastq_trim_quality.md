# fastq_trim_quality

Trim low-quality bases from FASTQ reads.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_trim_quality` |
| **Method ID** | `FastqParser.TrimByQuality` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Trims low-quality bases from both ends of FASTQ reads. The algorithm removes consecutive bases from the 5' and 3' ends that fall below the specified quality threshold.

Trimming process:
1. Scan from 5' end, removing bases until quality >= threshold
2. Scan from 3' end, removing bases until quality >= threshold
3. Return trimmed sequence with updated quality string

Reads that become empty after trimming are returned with empty sequence and quality.

## Core Documentation Reference

- Source: [FastqParser.cs#L264](../../../../Seqeron.Genomics/FastqParser.cs#L264)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | FASTQ format content to trim |
| `minQuality` | integer | No | `20` | Minimum quality score threshold |
| `encoding` | string | No | `"auto"` | Quality encoding: `"phred33"`, `"phred64"`, or `"auto"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries` | array | List of trimmed FASTQ records |
| `count` | integer | Number of records |
| `originalBases` | integer | Total bases before trimming |
| `trimmedBases` | integer | Total bases after trimming |
| `trimmedPercentage` | number | Percentage of bases removed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Minimum quality must be non-negative |

## Examples

### Example 1: Trim with default Q20 threshold

**Expected Tool Call:**
```json
{
  "tool": "fastq_trim_quality",
  "arguments": {
    "content": "@read1\nAAAATGCAAAA\n+\n!!!!!IIII!!"
  }
}
```

**Response:**
```json
{
  "entries": [
    {
      "id": "read1",
      "sequence": "TGC",
      "qualityString": "III",
      "length": 3
    }
  ],
  "count": 1,
  "originalBases": 11,
  "trimmedBases": 3,
  "trimmedPercentage": 72.7
}
```

### Example 2: Trim with custom threshold

**Expected Tool Call:**
```json
{
  "tool": "fastq_trim_quality",
  "arguments": {
    "content": "@read1\nATGCATGC\n+\nIIIIIIII",
    "minQuality": 30
  }
}
```

**Response:**
```json
{
  "entries": [
    {
      "id": "read1",
      "sequence": "ATGCATGC",
      "qualityString": "IIIIIIII",
      "length": 8
    }
  ],
  "count": 1,
  "originalBases": 8,
  "trimmedBases": 8,
  "trimmedPercentage": 0
}
```

## Performance

- **Time Complexity:** O(n) where n is total sequence length
- **Space Complexity:** O(n) for output records

## See Also

- [fastq_trim_adapter](fastq_trim_adapter.md) - Trim adapter sequences
- [fastq_filter](fastq_filter.md) - Filter by average quality
