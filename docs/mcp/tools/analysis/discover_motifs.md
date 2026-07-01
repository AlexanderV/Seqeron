# discover_motifs

De novo discovery of overrepresented k-mer motifs.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `discover_motifs` |
| **Method ID** | `MotifFinder.DiscoverMotifs` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Discovers overrepresented k-mers in a DNA sequence. Every length-k window is counted; k-mers
occurring at least `minCount` times are reported with their positions and an observed/expected
(O/E) enrichment. Under the i.i.d. uniform background, the expected count of a specific k-mer is
`E = (N − k + 1) / 4^k`, so `enrichment = observed / E` (Compeau & Pevzner). Values > 1 mean the
k-mer occurs more often than chance predicts. Input must be a valid DNA sequence.

## Core Documentation Reference

- Source: [MotifFinder.cs#L515](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs#L515)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `k` | integer | No | k-mer length (default 6, ≥ 1) |
| `minCount` | integer | No | Minimum occurrence count (default 2) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Motifs: `sequence, count, positions, enrichment` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Alternating repeat

**Input:** `{ "sequence": "ATATATAT", "k": 2, "minCount": 2 }`

N = 8, expected = 7/16 = 0.4375. AT occurs 4× (enrichment 9.143), TA occurs 3× (enrichment 6.857).

**Response (abridged):**
```json
{ "items": [
  { "sequence": "AT", "count": 4, "positions": [0,2,4,6], "enrichment": 9.1429 },
  { "sequence": "TA", "count": 3, "positions": [1,3,5], "enrichment": 6.8571 }
] }
```

### Example 2: High threshold filters everything

**Input:** `{ "sequence": "ACGTACGT", "k": 3, "minCount": 3 }`
→ no 3-mer occurs 3× → **Response:** `{ "items": [] }`

## Performance

- **Time Complexity:** O(n·k). **Space Complexity:** O(distinct k-mers).

## See Also

- [find_shared_motifs](find_shared_motifs.md)
- [most_frequent_kmers](most_frequent_kmers.md)
