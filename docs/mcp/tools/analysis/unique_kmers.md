# unique_kmers

The k-mers that occur exactly once (singletons).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `unique_kmers` |
| **Method ID** | `KmerAnalyzer.FindUniqueKmers` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **unique** (singleton) k-mers — those whose overlapping occurrence count
equals exactly 1. Note the distinction from *distinct*: a distinct-k-mer count counts
each different k-mer once regardless of multiplicity, whereas a *unique* k-mer must
appear only once. Counting is case-insensitive; order is unspecified. When `k`
exceeds the sequence length the result is empty.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L253](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L253)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `kmers` | array of string | k-mers occurring exactly once |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Singletons of ATGATG (k=3)

**User Prompt:**
> Which 3-mers appear only once in "ATGATG"?

**Expected Tool Call:**
```json
{
  "tool": "unique_kmers",
  "arguments": { "sequence": "ATGATG", "k": 3 }
}
```

**Response:**
```json
{ "kmers": ["TGA", "GAT"] }
```
ATG occurs twice (not unique); TGA and GAT once each.

### Example 2: No singletons (AAAA, k=2)

**User Prompt:**
> Unique 2-mers of "AAAA".

**Expected Tool Call:**
```json
{
  "tool": "unique_kmers",
  "arguments": { "sequence": "AAAA", "k": 2 }
}
```

**Response:**
```json
{ "kmers": [] }
```
AA appears 3 times, so there are no singletons.

## Performance

- **Time Complexity:** O(n) to build the k-mer table.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md) — full k-mer → count map
- [kmers_with_min_count](kmers_with_min_count.md) — count-thresholded k-mers
- [most_frequent_kmers](most_frequent_kmers.md) — the opposite extreme
