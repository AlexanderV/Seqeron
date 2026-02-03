# fastq_error_to_phred

Convert error probability to Phred quality score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_error_to_phred` |
| **Method ID** | `FastqParser.ErrorProbabilityToPhred` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Converts an error probability to its corresponding Phred quality score using the standard formula:

**Q = -10 * log10(P)**

Where:
- Q = Phred quality score (rounded to nearest integer)
- P = probability of incorrect base call

Special cases:
- P = 0 returns Q = 40 (maximum practical quality)
- P = 1 returns Q = 0

## Core Documentation Reference

- Source: [FastqParser.cs#L214](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/FastqParser.cs#L214)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `errorProbability` | number | Yes | Error probability (0-1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `errorProbability` | number | Input error probability |
| `phredScore` | integer | Calculated Phred quality score |

## Errors

| Code | Message |
|------|---------|
| 1001 | Error probability must be between 0 and 1 |

## Examples

### Example 1: Convert 1% error to Phred

**Expected Tool Call:**
```json
{
  "tool": "fastq_error_to_phred",
  "arguments": {
    "errorProbability": 0.01
  }
}
```

**Response:**
```json
{
  "errorProbability": 0.01,
  "phredScore": 20
}
```

### Example 2: Convert 0.1% error to Phred

**Expected Tool Call:**
```json
{
  "tool": "fastq_error_to_phred",
  "arguments": {
    "errorProbability": 0.001
  }
}
```

**Response:**
```json
{
  "errorProbability": 0.001,
  "phredScore": 30
}
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [fastq_phred_to_error](fastq_phred_to_error.md) - Convert Phred to error probability
- [fastq_statistics](fastq_statistics.md) - Calculate quality statistics
