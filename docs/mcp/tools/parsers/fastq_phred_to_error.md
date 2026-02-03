# fastq_phred_to_error

Convert Phred quality score to error probability.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Parsers |
| **Tool Name** | `fastq_phred_to_error` |
| **Method ID** | `FastqParser.PhredToErrorProbability` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Converts a Phred quality score to its corresponding error probability using the standard formula:

**P = 10^(-Q/10)**

Where:
- P = probability of incorrect base call
- Q = Phred quality score

Common conversions:
| Phred Score | Error Probability | Accuracy |
|-------------|-------------------|----------|
| 10 | 0.1 (1 in 10) | 90% |
| 20 | 0.01 (1 in 100) | 99% |
| 30 | 0.001 (1 in 1,000) | 99.9% |
| 40 | 0.0001 (1 in 10,000) | 99.99% |

## Core Documentation Reference

- Source: [FastqParser.cs#L206](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/FastqParser.cs#L206)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `phredScore` | integer | Yes | Phred quality score (non-negative) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `phredScore` | integer | Input Phred score |
| `errorProbability` | number | Probability of incorrect base call (0-1) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Phred score cannot be negative |

## Examples

### Example 1: Convert Q20 to probability

**Expected Tool Call:**
```json
{
  "tool": "fastq_phred_to_error",
  "arguments": {
    "phredScore": 20
  }
}
```

**Response:**
```json
{
  "phredScore": 20,
  "errorProbability": 0.01
}
```

### Example 2: Convert Q30 to probability

**Expected Tool Call:**
```json
{
  "tool": "fastq_phred_to_error",
  "arguments": {
    "phredScore": 30
  }
}
```

**Response:**
```json
{
  "phredScore": 30,
  "errorProbability": 0.001
}
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [fastq_error_to_phred](fastq_error_to_phred.md) - Convert error to Phred score
- [fastq_statistics](fastq_statistics.md) - Calculate quality statistics
