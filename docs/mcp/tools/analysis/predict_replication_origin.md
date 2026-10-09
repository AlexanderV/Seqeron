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
values, whether the signal is significant, and the full sets of minimizing/maximizing
prefix indices (Rosalind BA1F returns all minimizers, e.g. `53 97` for its sample). With
`circular: true` positions are reported modulo n.

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L247](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L247)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence, ideally a complete circular genome (min length 1) |
| `circular` | boolean | No | Treat the input as circular: prefix index n is the same junction as 0, positions reported mod n in [0, n−1] (default `false`: linear, [0, n]). Rotation-equivariant when total #G−#C = 0; otherwise the walk depends on the start |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `predictedOrigin` | integer | Prefix index of the cumulative-skew minimum |
| `predictedTerminus` | integer | Prefix index of the cumulative-skew maximum |
| `originSkew` | number | Cumulative skew at the origin |
| `terminusSkew` | number | Cumulative skew at the terminus |
| `isSignificant` | boolean | Threshold-free: true when max > min (non-zero amplitude) |
| `originPositions` | integer[] | All prefix indices minimizing the cumulative skew (Rosalind BA1F answer), ascending; first = `predictedOrigin` |
| `terminusPositions` | integer[] | All prefix indices maximizing the cumulative skew, ascending; first = `predictedTerminus` |

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
{ "predictedOrigin": 2, "predictedTerminus": 6, "originSkew": -2.0, "terminusSkew": 2.0, "isSignificant": true, "originPositions": [2], "terminusPositions": [6] }
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
{ "predictedOrigin": 0, "predictedTerminus": 3, "originSkew": 0.0, "terminusSkew": 3.0, "isSignificant": true, "originPositions": [0, 6], "terminusPositions": [3] }
```
With `"circular": true` the same input gives `"originPositions": [0]` (prefix 6 ≡ 0).

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(1) for the prediction plus O(#ties) for the position lists.

## See Also

- [cumulative_gc_skew](cumulative_gc_skew.md) — the underlying cumulative profile
- [gc_skew](gc_skew.md) — whole-sequence GC skew
