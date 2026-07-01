# count_kmers

Count every k-mer occurrence in a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `count_kmers` |
| **Method ID** | `KmerAnalyzer.CountKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Counts every overlapping k-mer (substring of length `k`) in a sequence and returns a
map of k-mer → occurrence count. Counting is case-insensitive. The counts sum to
`L − k + 1` (Wikipedia — K-mer). When `k` exceeds the sequence length the map is empty.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L20](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L20)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `counts` | object | Map of k-mer → occurrence count |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Homopolymer

**User Prompt:**
> Count 2-mers in "AAAA".

**Expected Tool Call:**
```json
{ "tool": "count_kmers", "arguments": { "sequence": "AAAA", "k": 2 } }
```

**Response:**
```json
{ "counts": { "AA": 3 } }
```
Overlapping windows AA, AA, AA ⇒ AA=3 (= L−k+1 = 3).

### Example 2: Mixed 3-mers

**Input:** `{ "sequence": "ATGATG", "k": 3 }`
→ **Response:** `{ "counts": { "ATG": 2, "TGA": 1, "GAT": 1 } }`

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(distinct k-mers).

## See Also

- [analyze_kmers](analyze_kmers.md)
- [kmer_frequencies](kmer_frequencies.md)
- [most_frequent_kmers](most_frequent_kmers.md)
