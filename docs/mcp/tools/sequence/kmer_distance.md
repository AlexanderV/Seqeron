# kmer_distance

Calculate k-mer based distance between two sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `kmer_distance` |
| **Method ID** | `KmerAnalyzer.KmerDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the k-mer distance between two sequences using Euclidean distance of k-mer frequencies. This alignment-free method is useful for comparing sequences of different lengths and provides a quick measure of sequence similarity. Lower values indicate more similar sequences; identical sequences have a distance of 0.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L165](../../../../src/Seqeron/Seqeron.Genomics/KmerAnalyzer.cs#L165)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First sequence (min length: 1) |
| `sequence2` | string | Yes | Second sequence (min length: 1) |
| `k` | integer | No | K-mer length (default: 3, minimum: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | number | K-mer distance (0 = identical k-mer composition) |
| `k` | integer | K-mer length used |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | K must be at least 1 |

## Examples

### Example 1: Compare identical sequences

**User Prompt:**
> What's the k-mer distance between "ATGCATGC" and itself?

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": {
    "sequence1": "ATGCATGC",
    "sequence2": "ATGCATGC",
    "k": 3
  }
}
```

**Response:**
```json
{
  "distance": 0,
  "k": 3
}
```

### Example 2: Compare different sequences

**User Prompt:**
> Compare "ATGCATGC" and "CCCCCCCC" using k-mer distance

**Expected Tool Call:**
```json
{
  "tool": "kmer_distance",
  "arguments": {
    "sequence1": "ATGCATGC",
    "sequence2": "CCCCCCCC",
    "k": 3
  }
}
```

**Response:**
```json
{
  "distance": 0.87,
  "k": 3
}
```

## Performance

- **Time Complexity:** O(n + m) where n, m are sequence lengths
- **Space Complexity:** O(4^k) for k-mer frequency storage

## See Also

- [kmer_count](kmer_count.md) - Count k-mer frequencies
- [kmer_analyze](kmer_analyze.md) - Comprehensive k-mer statistics
- [edit_distance](../core/edit_distance.md) - Levenshtein distance
