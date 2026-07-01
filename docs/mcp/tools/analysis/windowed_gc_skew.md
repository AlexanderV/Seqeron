# windowed_gc_skew

Sliding-window GC skew along a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `windowed_gc_skew` |
| **Method ID** | `GcSkewCalculator.CalculateWindowedGcSkew` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **GC skew in sliding windows** along a DNA sequence. Each point reports
the window skew `(G − C) / (G + C)`, the window center position (`i + windowSize/2`),
and the window bounds `[i, i + windowSize − 1]`. Windows are advanced by `stepSize`.

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L73](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L73)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `windowSize` | integer | No | Window size in bp (default 1000, ≥ 1) |
| `stepSize` | integer | No | Step size in bp (default 100, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ position, gcSkew, windowStart, windowEnd }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |
| 1002 | Step size must be at least 1 |

## Examples

### Example 1: Two non-overlapping windows

**User Prompt:**
> Windowed GC skew of "GGGGCCCC" with window 4, step 4.

**Expected Tool Call:**
```json
{
  "tool": "windowed_gc_skew",
  "arguments": { "sequence": "GGGGCCCC", "windowSize": 4, "stepSize": 4 }
}
```

**Response:**
```json
{ "items": [ { "position": 2, "gcSkew": 1.0, "windowStart": 0, "windowEnd": 3 }, { "position": 6, "gcSkew": -1.0, "windowStart": 4, "windowEnd": 7 } ] }
```

### Example 2: Window larger than sequence

**User Prompt:**
> Windowed GC skew of "GC" with window 4.

**Expected Tool Call:**
```json
{
  "tool": "windowed_gc_skew",
  "arguments": { "sequence": "GC", "windowSize": 4, "stepSize": 4 }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n · windowSize / stepSize).
- **Space Complexity:** O(number of windows).

## See Also

- [gc_skew](gc_skew.md) — whole-sequence GC skew
- [cumulative_gc_skew](cumulative_gc_skew.md) — cumulative skew
- [predict_replication_origin](predict_replication_origin.md) — origin/terminus from skew
