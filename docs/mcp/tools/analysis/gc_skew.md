# gc_skew

Whole-sequence GC skew.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `gc_skew` |
| **Method ID** | `GcSkewCalculator.CalculateGcSkew` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **GC skew** of a DNA sequence: `(G − C) / (G + C)`. The value lies in
`[−1, 1]`; it is 0 when the sequence contains no G or C. GC skew shifts sign around the
replication origin and terminus in bacterial genomes.

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L30](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L30)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `gcSkew` | number | GC skew in `[−1, 1]` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: All G (skew = 1)

**User Prompt:**
> GC skew of "GGGG".

**Expected Tool Call:**
```json
{
  "tool": "gc_skew",
  "arguments": { "sequence": "GGGG" }
}
```

**Response:**
```json
{ "gcSkew": 1.0 }
```

### Example 2: Balanced (skew = 0)

**User Prompt:**
> GC skew of "GGGGCCCC".

**Expected Tool Call:**
```json
{
  "tool": "gc_skew",
  "arguments": { "sequence": "GGGGCCCC" }
}
```

**Response:**
```json
{ "gcSkew": 0.0 }
```

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(1).

## See Also

- [windowed_gc_skew](windowed_gc_skew.md) — sliding-window GC skew
- [cumulative_gc_skew](cumulative_gc_skew.md) — cumulative skew for origin finding
- [at_skew](at_skew.md) — AT skew
