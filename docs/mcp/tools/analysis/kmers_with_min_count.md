# kmers_with_min_count

k-mers occurring at least `minCount` times, ordered by count descending.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `kmers_with_min_count` |
| **Method ID** | `KmerAnalyzer.FindKmersWithMinCount` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **recurrent** k-mers — those whose overlapping occurrence count is at
least `minCount` (`Count(Text, Pattern) ≥ t`, per Compeau & Pevzner) — as
`(kmer, count)` pairs ordered by count descending. With `minCount ≤ 1` every distinct
k-mer qualifies. Counting is case-insensitive; when `k` exceeds the sequence length
the result is empty.

## Core Documentation Reference

- Source: [KmerAnalyzer.cs#L274](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs#L274)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `k` | integer | Yes | k-mer length (> 0) |
| `minCount` | integer | Yes | Inclusive minimum occurrence count |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `(kmer, count)` pairs with count ≥ minCount, ordered by count descending |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | k must be positive |

## Examples

### Example 1: Recurrent trimer (minCount = 2)

**User Prompt:**
> Which 3-mers occur at least twice in "ATGATG"?

**Expected Tool Call:**
```json
{
  "tool": "kmers_with_min_count",
  "arguments": { "sequence": "ATGATG", "k": 3, "minCount": 2 }
}
```

**Response:**
```json
{ "items": [ { "kmer": "ATG", "count": 2 } ] }
```

### Example 2: All distinct trimers (minCount = 1)

**User Prompt:**
> List every 3-mer of "ATGATG" with its count.

**Expected Tool Call:**
```json
{
  "tool": "kmers_with_min_count",
  "arguments": { "sequence": "ATGATG", "k": 3, "minCount": 1 }
}
```

**Response:**
```json
{ "items": [ { "kmer": "ATG", "count": 2 }, { "kmer": "TGA", "count": 1 }, { "kmer": "GAT", "count": 1 } ] }
```
ATG (count 2) sorts before the count-1 k-mers.

## Performance

- **Time Complexity:** O(n) to build the k-mer table + O(d log d) to sort `d` distinct k-mers.
- **Space Complexity:** O(distinct k-mers).

## See Also

- [count_kmers](count_kmers.md) — full k-mer → count map
- [most_frequent_kmers](most_frequent_kmers.md) — only the top-count k-mers
- [unique_kmers](unique_kmers.md) — the singletons
