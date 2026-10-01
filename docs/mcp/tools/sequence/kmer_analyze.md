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

- Source: [KmerAnalyzer.cs#L654](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L654)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |
| `k` | integer | No | K-mer length (default: 3, minimum: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `totalKmers` | integer | Total number of k-mers including multiplicity (Jellyfish Total), `L − k + 1` |
| `uniqueKmers` | integer | Number of **distinct** k-mers (legacy name; = `distinctKmers`; not Jellyfish "Unique") |
| `distinctKmers` | integer | Number of distinct k-mers (Jellyfish Distinct) |
| `singletonKmers` | integer | Number of k-mers occurring exactly once (Jellyfish Unique) |
| `maxCount` | integer | Maximum frequency of any k-mer |
| `minCount` | integer | Minimum frequency of any k-mer |
| `averageCount` | number | Exact mean multiplicity `total/distinct` |
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
  "entropy": 1.970950594454669,
  "k": 3,
  "distinctKmers": 4,
  "singletonKmers": 0
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
