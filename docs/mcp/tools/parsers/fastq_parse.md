# fastq_parse

Parse FASTQ format content into sequence entries with quality scores.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_parse` |
| **Method ID** | `FastqParser.Parse` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Parses FASTQ format string into individual sequence entries with quality scores. Each entry includes the sequence ID, optional description, DNA sequence, quality string, and decoded quality scores. FASTQ is the standard format for sequencing reads with quality information.

Supports multiple quality encodings:
- **Phred+33** (Sanger/Illumina 1.8+): ASCII offset 33, scores 0-93
- **Phred+64** (Illumina 1.3-1.7): ASCII offset 64, scores 0-62
- **Auto**: Automatically detects encoding based on quality characters

## Core Documentation Reference

- Source: [FastqParser.cs#L21](../../../../src/Seqeron/Seqeron.Genomics/FastqParser.cs#L21)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | FASTQ format content to parse |
| `encoding` | string | No | `"auto"` | Quality encoding: `"phred33"`, `"phred64"`, or `"auto"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `entries` | array | List of parsed FASTQ records |
| `entries[].id` | string | Sequence identifier (first word of header) |
| `entries[].description` | string? | Optional description (rest of header) |
| `entries[].sequence` | string | DNA sequence |
| `entries[].qualityString` | string | Raw quality string |
| `entries[].qualityScores` | integer[] | Decoded Phred quality scores |
| `entries[].length` | integer | Sequence length in base pairs |
| `count` | integer | Number of records parsed |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto' |

## Examples

### Example 1: Parse FASTQ with auto encoding

**User Prompt:**
> Parse this FASTQ content with two reads

**Expected Tool Call:**
```json
{
  "tool": "fastq_parse",
  "arguments": {
    "content": "@read1 sample1\nATGCATGC\n+\nIIIIIIII\n@read2\nGGGCCC\n+\nHHHHHH"
  }
}
```

**Response:**
```json
{
  "entries": [
    {
      "id": "read1",
      "description": "sample1",
      "sequence": "ATGCATGC",
      "qualityString": "IIIIIIII",
      "qualityScores": [40, 40, 40, 40, 40, 40, 40, 40],
      "length": 8
    },
    {
      "id": "read2",
      "description": null,
      "sequence": "GGGCCC",
      "qualityString": "HHHHHH",
      "qualityScores": [39, 39, 39, 39, 39, 39],
      "length": 6
    }
  ],
  "count": 2
}
```

### Example 2: Parse with explicit Phred+33 encoding

**User Prompt:**
> Parse this Illumina FASTQ data using Phred+33 encoding

**Expected Tool Call:**
```json
{
  "tool": "fastq_parse",
  "arguments": {
    "content": "@SRR001\nACGT\n+\n!!!!",
    "encoding": "phred33"
  }
}
```

**Response:**
```json
{
  "entries": [
    {
      "id": "SRR001",
      "description": null,
      "sequence": "ACGT",
      "qualityString": "!!!!",
      "qualityScores": [0, 0, 0, 0],
      "length": 4
    }
  ],
  "count": 1
}
```

## Performance

- **Time Complexity:** O(n) where n is total content length
- **Space Complexity:** O(n) for storing parsed records

## See Also

- [fasta_parse](fasta_parse.md) - Parse FASTA format (without quality)
- [fastq_statistics](fastq_statistics.md) - Calculate quality statistics
- [fastq_filter](fastq_filter.md) - Filter reads by quality
