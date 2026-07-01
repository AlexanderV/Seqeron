# entropy_profile

Shannon entropy in sliding windows along a sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `entropy_profile` |
| **Method ID** | `SequenceStatistics.CalculateEntropyProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Produces the Shannon-entropy profile of a sequence: for each sliding window of length
`windowSize` (advanced by `stepSize`) it reports `H = −Σ pᵢ log₂ pᵢ` in bits over the window's
letter frequencies (case-insensitive; non-letters ignored). Windows are emitted for offsets
`0, stepSize, 2·stepSize, …` up to `length − windowSize`. Maximum entropy is `log₂ k` for k
distinct symbols (2 bits for the 4-letter DNA alphabet).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L963](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L963)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Sequence to analyze (min length 1) |
| `windowSize` | integer | No | Window size (default 50, ≥ 1) |
| `stepSize` | integer | No | Step size (default 1, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `values` | number[] | Per-window entropy in bits |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Window size must be at least 1 |
| 1003 | Step size must be at least 1 |

## Examples

### Example 1: Homopolymer (zero entropy)

**Input:** `{ "sequence": "AAAA", "windowSize": 2, "stepSize": 1 }`
→ three "AA" windows → **Response:** `{ "values": [0, 0, 0] }`

### Example 2: Uniform window (max entropy)

**Input:** `{ "sequence": "ATGC", "windowSize": 4 }`
→ one window with equal base frequencies → **Response:** `{ "values": [2.0] }`

## Performance

- **Time Complexity:** O(n · windowSize). **Space Complexity:** O(n / stepSize).

## See Also

- [windowed_complexity](windowed_complexity.md)
- [gc_content_profile](gc_content_profile.md)
