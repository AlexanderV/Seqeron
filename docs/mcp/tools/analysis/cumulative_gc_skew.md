# cumulative_gc_skew

Cumulative GC skew along a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `cumulative_gc_skew` |
| **Method ID** | `GcSkewCalculator.CalculateCumulativeGcSkew` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the cumulative GC skew across a DNA sequence using non-overlapping windows
(step size equals the window size). For each window the per-window skew is `(G − C)/(G + C)`
and the cumulative value is the running sum of those window skews. The reported position is
the window center `windowStart + windowSize/2`. On a complete circular bacterial genome the
**minimum** of the cumulative curve approximates the replication origin and the **maximum**
the terminus (Grigoriev 1998; Lobry 1996).

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L115](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L115)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `windowSize` | integer | No | Window size (default 1000, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | One point per window: `position`, `gcSkew`, `cumulativeGcSkew` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1003 | Window size must be at least 1 |

## Examples

### Example 1: G-block then C-block

**Input:** `{ "sequence": "GGGGCCCC", "windowSize": 4 }`

Window 0 "GGGG" skew +1 (cumulative 1); window 1 "CCCC" skew −1 (cumulative 0).

**Response:**
```json
{ "items": [
  { "position": 2, "gcSkew": 1.0, "cumulativeGcSkew": 1.0 },
  { "position": 6, "gcSkew": -1.0, "cumulativeGcSkew": 0.0 }
] }
```

### Example 2: Short sequence, no full window

**Input:** `{ "sequence": "GC", "windowSize": 4 }`
→ **Response:** `{ "items": [] }` (no complete window fits).

## Performance

- **Time Complexity:** O(n). **Space Complexity:** O(n / windowSize).

## See Also

- [gc_skew](gc_skew.md)
- [windowed_gc_skew](windowed_gc_skew.md)
- [predict_replication_origin](predict_replication_origin.md)
