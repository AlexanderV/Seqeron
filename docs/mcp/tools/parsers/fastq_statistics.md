# fastq_statistics

Calculate quality statistics for FASTQ data.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_statistics` |
| **Method ID** | `FastqParser.CalculateStatistics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates comprehensive quality statistics for FASTQ sequencing data. Returns read counts, base counts, length statistics, quality metrics (Q20/Q30 percentages), and GC content. Essential for quality control of sequencing runs.

## Core Documentation Reference

- Source: [FastqParser.cs#L344](../../../../src/Seqeron/Seqeron.Genomics/FastqParser.cs#L344)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `content` | string | Yes | - | FASTQ format content to analyze |
| `encoding` | string | No | `"auto"` | Quality encoding: `"phred33"`, `"phred64"`, or `"auto"` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalReads` | integer | Total number of reads |
| `totalBases` | integer | Total number of bases across all reads |
| `meanReadLength` | number | Average read length |
| `meanQuality` | number | Mean Phred quality score |
| `minReadLength` | integer | Minimum read length |
| `maxReadLength` | integer | Maximum read length |
| `q20Percentage` | number | Percentage of bases with Q >= 20 |
| `q30Percentage` | number | Percentage of bases with Q >= 30 |
| `gcContent` | number | GC content as fraction (0.0-1.0) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Content cannot be null or empty |
| 1002 | Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto' |

## Examples

### Example 1: Basic quality statistics

**User Prompt:**
> What are the quality statistics for this FASTQ data?

**Expected Tool Call:**
```json
{
  "tool": "fastq_statistics",
  "arguments": {
    "content": "@read1\nATGCATGC\n+\nIIIIIIII\n@read2\nGGGCCC\n+\nHHHHHH"
  }
}
```

**Response:**
```json
{
  "totalReads": 2,
  "totalBases": 14,
  "meanReadLength": 7.0,
  "meanQuality": 39.57,
  "minReadLength": 6,
  "maxReadLength": 8,
  "q20Percentage": 100.0,
  "q30Percentage": 100.0,
  "gcContent": 0.5
}
```

### Example 2: Statistics with explicit encoding

**User Prompt:**
> Calculate statistics for this Illumina 1.8+ FASTQ using Phred+33 encoding

**Expected Tool Call:**
```json
{
  "tool": "fastq_statistics",
  "arguments": {
    "content": "@SRR001\nACGTACGT\n+\n!!!!!!!!",
    "encoding": "phred33"
  }
}
```

**Response:**
```json
{
  "totalReads": 1,
  "totalBases": 8,
  "meanReadLength": 8.0,
  "meanQuality": 0.0,
  "minReadLength": 8,
  "maxReadLength": 8,
  "q20Percentage": 0.0,
  "q30Percentage": 0.0,
  "gcContent": 0.5
}
```

## Performance

- **Time Complexity:** O(n) where n is total content length
- **Space Complexity:** O(n) for storing parsed records during calculation

## Quality Metrics

- **Q20**: 1% error rate (99% accuracy)
- **Q30**: 0.1% error rate (99.9% accuracy)
- **GC Content**: Fraction of G and C bases, useful for detecting contamination

## See Also

- [fastq_parse](fastq_parse.md) - Parse FASTQ into records
- [fastq_filter](fastq_filter.md) - Filter reads by quality
