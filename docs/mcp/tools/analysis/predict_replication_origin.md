# predict_replication_origin

Predict replication origin and terminus from cumulative GC skew.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_replication_origin` |
| **Method ID** | `GcSkewCalculator.PredictReplicationOrigin` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts the **replication origin and terminus** from cumulative GC-skew extrema: the
origin is approximated by the prefix index where the cumulative skew is minimal, and
the terminus by the index where it is maximal (Lobry 1996; Grigoriev 1998). Works best
on complete circular bacterial genomes. Returns both predicted positions, their skew
values, and whether the signal is significant.

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L247](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L247)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence, ideally a complete circular genome (min length 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `predictedOrigin` | integer | Prefix index of the cumulative-skew minimum |
| `predictedTerminus` | integer | Prefix index of the cumulative-skew maximum |
| `originSkew` | number | Cumulative skew at the origin |
| `terminusSkew` | number | Cumulative skew at the terminus |
| `isSignificant` | boolean | Whether the skew signal is significant |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: CCGGGG

**User Prompt:**
> Predict the replication origin of "CCGGGG".

**Expected Tool Call:**
```json
{
  "tool": "predict_replication_origin",
  "arguments": { "sequence": "CCGGGG" }
}
```

**Response:**
```json
{ "predictedOrigin": 2, "predictedTerminus": 6, "originSkew": -2.0, "terminusSkew": 2.0 }
```
The cumulative skew reaches its minimum (−2) after the two leading C's and its maximum
(+2) at the end.

### Example 2: GGGCCC

**User Prompt:**
> Predict the replication origin of "GGGCCC".

**Expected Tool Call:**
```json
{
  "tool": "predict_replication_origin",
  "arguments": { "sequence": "GGGCCC" }
}
```

**Response:**
```json
{ "predictedOrigin": 0, "predictedTerminus": 3, "originSkew": 0.0, "terminusSkew": 3.0 }
```

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(n) for the cumulative profile.

## See Also

- [cumulative_gc_skew](cumulative_gc_skew.md) — the underlying cumulative profile
- [gc_skew](gc_skew.md) — whole-sequence GC skew
