# fastq_trim_adapter

Trim adapter sequences from FASTQ reads.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_trim_adapter` |
| **Method ID** | `FastqParser.TrimAdapter` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Removes adapter sequences from FASTQ reads. The algorithm searches for the adapter at the 3' end of the sequence and within the sequence, trimming everything from the adapter start position.

Search strategy:
1. Check for full or partial adapter match at 3' end (overlap >= minOverlap)
2. Check for full adapter match anywhere in the sequence
3. If found, trim the read at the adapter start position

Case-insensitive matching is used for adapter detection.

## Core Documentation Reference

- Source: [FastqParser.cs#L297](../../../../Seqeron.Genomics/FastqParser.cs#L297)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | FASTQ format content to trim |
| `adapter` | string | Yes | - | Adapter sequence to remove |
| `minOverlap` | integer | No | `5` | Minimum overlap length to consider a match |
| `encoding` | string | No | `"auto"` | Quality encoding: `"phred33"`, `"phred64"`, or `"auto"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries` | array | List of trimmed FASTQ records |
| `count` | integer | Number of records |
| `readsWithAdapter` | integer | Number of reads with adapter found |
| `originalBases` | integer | Total bases before trimming |
| `trimmedBases` | integer | Total bases after trimming |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Adapter cannot be null or empty |
| 1003 | Minimum overlap must be at least 1 |

## Examples

### Example 1: Trim Illumina adapter

**Expected Tool Call:**
```json
{
  "tool": "fastq_trim_adapter",
  "arguments": {
    "content": "@read1\nATGCATGCAAAAAAAA\n+\nIIIIIIIIIIIIIIII",
    "adapter": "AAAAAAAA"
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
  "readsWithAdapter": 1,
  "originalBases": 16,
  "trimmedBases": 8
}
```

### Example 2: No adapter found

**Expected Tool Call:**
```json
{
  "tool": "fastq_trim_adapter",
  "arguments": {
    "content": "@read1\nATGCATGC\n+\nIIIIIIII",
    "adapter": "GGGGGGGG"
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
  "readsWithAdapter": 0,
  "originalBases": 8,
  "trimmedBases": 8
}
```

## Performance

- **Time Complexity:** O(n * m) where n is sequence length and m is adapter length
- **Space Complexity:** O(n) for output records

## See Also

- [fastq_trim_quality](fastq_trim_quality.md) - Trim by quality score
- [fastq_filter](fastq_filter.md) - Filter by average quality
