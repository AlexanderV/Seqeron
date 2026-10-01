# windowed_complexity

Sliding-window Shannon entropy and linguistic complexity (DNA).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `windowed_complexity` |
| **Method ID** | `SequenceComplexity.CalculateWindowedComplexity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes, in sliding windows along a DNA sequence, both the **Shannon entropy**
(bits/symbol) and the **linguistic complexity**. Each point reports the window center
position (`start + windowSize/2`), the inclusive window bounds `[start, start+w-1]`,
the entropy and the linguistic complexity (word lengths 1..min(`lcMaxWordLength`, w), default 6). Windows advance by `stepSize`.
N / IUPAC codes are accepted; a window containing any non-ACGT symbol is skipped (no point), as BBDuk scores only windows with `ns() < 1`.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L470](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L470)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence, A/C/G/T + IUPAC codes (min length 1) |
| `windowSize` | integer | No | Window size (default 64) |
| `stepSize` | integer | No | Step size (default 10) |
| `lcMaxWordLength` | integer | No | Per-window LC word-length cap (default 6; ≥ windowSize gives the all-length LC, e.g. `ACGTACGT` → 13/16) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ position, shannonEntropy, linguisticComplexity, windowStart, windowEnd }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Three non-overlapping windows

**User Prompt:**
> Windowed complexity of "ACGTACGTAAAAAAAAACGTACGT" with window 8, step 8.

**Expected Tool Call:**
```json
{
  "tool": "windowed_complexity",
  "arguments": { "sequence": "ACGTACGTAAAAAAAAACGTACGT", "windowSize": 8, "stepSize": 8 }
}
```

**Response:**
```json
{ "items": [ { "position": 4, "shannonEntropy": 2.0, "linguisticComplexity": 0.7931034482758621, "windowStart": 0, "windowEnd": 7 }, { "position": 12, "shannonEntropy": 0.0, "linguisticComplexity": 0.20689655172413793, "windowStart": 8, "windowEnd": 15 }, { "position": 20, "shannonEntropy": 2.0, "linguisticComplexity": 0.7931034482758621, "windowStart": 16, "windowEnd": 23 } ] }
```
The uniform ACGTACGT windows have maximal entropy (log₂4 = 2.0); the homopolymer
window has entropy 0. Linguistic complexity (summation form, word lengths 1..min(6, w)):
ACGTACGT = 23/29 = 0.7931034482758621, AAAAAAAA = 6/29 = 0.20689655172413793.

## Performance

- **Time Complexity:** O(n · windowSize / stepSize).
- **Space Complexity:** O(number of windows).

## See Also

- [find_low_complexity_regions](find_low_complexity_regions.md)
- [entropy_profile](entropy_profile.md)
- [dust_score](dust_score.md)
