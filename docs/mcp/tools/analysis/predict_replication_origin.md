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

With `windowSize` (Core `PredictReplicationOrigin(DnaSequence, int windowSize)`) the prediction uses
Grigoriev's (1998) **windowed** cumulative skew diagram instead: the running sum of `(G−C)/(G+C)` over
adjacent, non-overlapping complete windows; origin/terminus are the centres (`start + windowSize/2`) of
the first windows with the minimum/maximum running sum, and `originSkew`/`terminusSkew` are those sums
(equal to Biopython `numpy.cumsum(GC_skew(seq, w)[:n//w])` argmin/argmax). The windowed prediction is
linear only: `circular: true` together with `windowSize` is rejected (error 1001).
`originPositions`/`terminusPositions` are per-base prefix indices and are empty in windowed mode.

With `skewIndexWindow` (Core `CalculateSkewIndex(DnaSequence, int windowSize)`) the response also
carries `skewIndex`, the SkewIT Skew Index of Lu & Salzberg (2020, PLoS Comput Biol 16:e1008439),
computed exactly as the authors' `skewi.py` (SkewIT default window 20000). It is `null` when skewi.py
prints no value (fewer than 13 windows, or no window with #G ≠ #C) and when `skewIndexWindow` is
omitted. There is no universal SkewI cutoff; SkewIT publishes per-genus thresholds (e.g. Escherichia
0.7110) and `isSignificant` is unrelated to it. skewi.py's input filters (≥ 500 kb, "complete" header,
no plasmids) are not applied.

## Core Documentation Reference

- Source: [GcSkewCalculator.cs#L657](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L657)
- Windowed (`windowSize`): [GcSkewCalculator.cs#L693](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L693); Skew Index (`skewIndexWindow`): [GcSkewCalculator.cs#L769](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs#L769)
- Evidence: `docs/Evidence/SEQ-REPLICATION-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence, ideally a complete circular genome (min length 1) |
| `circular` | boolean | No | Treat the input as circular: prefix index n is the same junction as 0, positions reported mod n in [0, n−1] (default `false`: linear, [0, n]). Rotation-equivariant when total #G−#C = 0; otherwise the walk depends on the start. Not combinable with `windowSize` |
| `windowSize` | integer | No | Grigoriev window in bp (≥ 1). Omit (default) for the per-base walk |
| `skewIndexWindow` | integer | No | SkewIT window k in bp (≥ 1; SkewIT default 20000). Adds `skewIndex`. Omit to skip |

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
| `skewIndex` | number \| null | SkewIT Skew Index in (0, 1] (only with `skewIndexWindow`; null when skewi.py reports none) |

In windowed mode `predictedOrigin`/`predictedTerminus` are window centres, the skews are cumulative
window sums, and both position arrays are empty.

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1001 | circular is not supported with windowSize (the windowed Grigoriev prediction is linear) |
| 1001 | Window size must be at least 1 / Skew index window must be at least 1 |

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

### Example 3: Windowed prediction and Skew Index (Rosalind BA1F sample)

**Expected Tool Call:**
```json
{
  "tool": "predict_replication_origin",
  "arguments": { "sequence": "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG", "windowSize": 10, "skewIndexWindow": 4 }
}
```

**Response:**
```json
{ "predictedOrigin": 45, "predictedTerminus": 15, "originSkew": -0.3095238095238095, "terminusSkew": 0.5, "isSignificant": true, "originPositions": [], "terminusPositions": [], "skewIndex": 0.16 }
```
Biopython 1.88 `numpy.cumsum(GC_skew(seq, 10)[:10])`: minimum −0.3095… at window 4 (centre 45),
maximum 0.5 at window 1 (centre 15). SkewIT `skewi.py -k 4 --min-len 0` prints 0.16.

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(1) for the prediction plus O(#ties) for the position lists.
- **Skew Index:** O(L · 0.08 L) for L = ⌈n/k⌉ windows (skewi.py's split scan, via prefix sums).

## See Also

- [cumulative_gc_skew](cumulative_gc_skew.md) — the underlying cumulative profile
- [gc_skew](gc_skew.md) — whole-sequence GC skew
