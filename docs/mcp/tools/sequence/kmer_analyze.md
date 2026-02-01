# kmer_analyze

Comprehensive k-mer analysis of a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `kmer_analyze` |
| **Method ID** | `KmerAnalyzer.AnalyzeKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs comprehensive k-mer analysis on a sequence, returning statistics about the frequency distribution including total count, unique count, min/max/average frequencies, and Shannon entropy. This provides a complete picture of the k-mer composition of a sequence.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L363](../../../../src/Seqeron/Seqeron.Genomics/KmerAnalyzer.cs#L363)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |
| `k` | integer | No | K-mer length (default: 3, minimum: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalKmers` | integer | Total number of k-mers in sequence |
| `uniqueKmers` | integer | Number of distinct k-mers |
| `maxCount` | integer | Maximum frequency of any k-mer |
| `minCount` | integer | Minimum frequency of any k-mer |
| `averageCount` | number | Average k-mer frequency |
| `entropy` | number | Shannon entropy of k-mer distribution (bits) |
| `k` | integer | K-mer length used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | K must be at least 1 |

## Examples

### Example 1: Analyze DNA sequence

**User Prompt:**
> Analyze the 3-mer composition of "ATGCATGCATGC"

**Expected Tool Call:**
```json
{
  "tool": "kmer_analyze",
  "arguments": {
    "sequence": "ATGCATGCATGC",
    "k": 3
  }
}
```

**Response:**
```json
{
  "totalKmers": 10,
  "uniqueKmers": 4,
  "maxCount": 3,
  "minCount": 2,
  "averageCount": 2.5,
  "entropy": 1.97,
  "k": 3
}
```

### Example 2: Analyze with different k

**User Prompt:**
> What's the 4-mer statistics for "ATGCATGCATGC"?

**Expected Tool Call:**
```json
{
  "tool": "kmer_analyze",
  "arguments": {
    "sequence": "ATGCATGCATGC",
    "k": 4
  }
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(4^k) for storing unique k-mers

## See Also

- [kmer_count](kmer_count.md) - Get individual k-mer counts
- [kmer_entropy](kmer_entropy.md) - K-mer entropy only
- [kmer_distance](kmer_distance.md) - Compare sequences by k-mers
