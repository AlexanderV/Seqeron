# gc_content_profile

GC content in sliding windows along a nucleotide sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `gc_content_profile` |
| **Method ID** | `SequenceStatistics.CalculateGcContentProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **GC-content profile**: the GC percentage of each sliding window of
`windowSize` bases, advanced by `stepSize`. GC content is
`(G + C) / (A + T + U + G + C) × 100` over the standard bases in the window (ambiguous
symbols such as N are excluded, matching Biopython `gc_fraction`); U counts as a
non-GC base. When the window exceeds the sequence length the result is empty.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L905](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L905)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence (min length 1) |
| `windowSize` | integer | No | Window size (default 100, ≥ 1) |
| `stepSize` | integer | No | Step size (default 1, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `values` | array of number | GC percentage per window (0–100) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |
| 1002 | Step size must be at least 1 |

## Examples

### Example 1: GC-rich windows (window 2, step 1)

**User Prompt:**
> GC profile of "GGCC" with window 2.

**Expected Tool Call:**
```json
{
  "tool": "gc_content_profile",
  "arguments": { "sequence": "GGCC", "windowSize": 2, "stepSize": 1 }
}
```

**Response:**
```json
{ "values": [100.0, 100.0, 100.0] }
```
Windows GG, GC, CC are all 100% GC.

### Example 2: AT-only windows, non-overlapping (window 2, step 2)

**User Prompt:**
> GC profile of "ATAT" with window 2, step 2.

**Expected Tool Call:**
```json
{
  "tool": "gc_content_profile",
  "arguments": { "sequence": "ATAT", "windowSize": 2, "stepSize": 2 }
}
```

**Response:**
```json
{ "values": [0.0, 0.0] }
```
Both non-overlapping windows (AT, AT) contain no G or C.

## Performance

- **Time Complexity:** O(n · windowSize / stepSize).
- **Space Complexity:** O(number of windows).

## See Also

- [gc_content](../sequence/gc_content.md) — whole-sequence GC content
- [analyze_gc_content](analyze_gc_content.md) — comprehensive GC report
- [entropy_profile](entropy_profile.md) — windowed Shannon entropy
