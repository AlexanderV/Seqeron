# find_clumps

Find k-mers that clump within a sliding window.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_clumps` |
| **Method ID** | `KmerAnalyzer.FindClumps` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds all k-mers that occur at least `minOccurrences` times within **some** sliding window of
size `windowSize` (the Clump Finding Problem, Compeau & Pevzner Ch. 1). Counting is
case-insensitive and overlapping. The returned set is unordered. `windowSize` must be at least
`k`.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L356](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L356)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to scan (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `windowSize` | integer | Yes | Sliding window size (≥ k) |
| `minOccurrences` | integer | Yes | Minimum occurrences within a window (> 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `kmers` | string[] | k-mers forming a clump (unordered) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |
| 1003 | Window size must be at least k |
| 1003 | Minimum occurrences must be positive |

## Examples

### Example 1: Homopolymer clump

**Input:** `{ "sequence": "AAAA", "k": 2, "windowSize": 4, "minOccurrences": 3 }`
→ window "AAAA" has AA ×3 → **Response:** `{ "kmers": ["AA"] }`

### Example 2: Two clumps

**Input:** `{ "sequence": "AAAACCCC", "k": 2, "windowSize": 4, "minOccurrences": 3 }`
→ **Response:** `{ "kmers": ["AA", "CC"] }` (order unspecified).

## Performance

- **Time Complexity:** O(n · windowSize). **Space Complexity:** O(distinct k-mers per window).

## See Also

- [most_frequent_kmers](most_frequent_kmers.md)
- [kmer_positions](kmer_positions.md)
