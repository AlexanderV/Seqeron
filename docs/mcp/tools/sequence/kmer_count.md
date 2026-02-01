# kmer_count

Count k-mer frequencies in a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `kmer_count` |
| **Method ID** | `KmerAnalyzer.CountKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts all k-mers (substrings of length k) in a sequence and returns their frequencies. K-mer analysis is fundamental in bioinformatics for tasks like genome assembly, sequence comparison, and motif discovery.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L20](../../../../src/Seqeron/Seqeron.Genomics/KmerAnalyzer.cs#L20)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to analyze (min length: 1) |
| `k` | integer | No | K-mer length (default: 3, minimum: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `counts` | object | Dictionary mapping k-mers to their counts |
| `k` | integer | K-mer length used |
| `uniqueKmers` | integer | Number of unique k-mers found |
| `totalKmers` | integer | Total k-mer count |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | K must be at least 1 |

## Examples

### Example 1: Count 3-mers in DNA

**User Prompt:**
> Count all 3-mers in "ATGCATGC"

**Expected Tool Call:**
```json
{
  "tool": "kmer_count",
  "arguments": {
    "sequence": "ATGCATGC",
    "k": 3
  }
}
```

**Response:**
```json
{
  "counts": {
    "ATG": 2,
    "TGC": 2,
    "GCA": 1,
    "CAT": 1
  },
  "k": 3,
  "uniqueKmers": 4,
  "totalKmers": 6
}
```

### Example 2: Count dinucleotides

**User Prompt:**
> What dinucleotides are in "ATGCATGC"?

**Expected Tool Call:**
```json
{
  "tool": "kmer_count",
  "arguments": {
    "sequence": "ATGCATGC",
    "k": 2
  }
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(4^k) for storing unique k-mers

## See Also

- [kmer_analyze](kmer_analyze.md) - Comprehensive k-mer statistics
- [kmer_entropy](kmer_entropy.md) - K-mer entropy calculation
- [kmer_distance](kmer_distance.md) - K-mer based distance
