# reversal_distance

Lower-bound reversal distance between two permutations.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `reversal_distance` |
| **Method ID** | `ComparativeGenomics.CalculateReversalDistance` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes a **lower bound on the reversal (inversion) distance** between two
equal-length gene-order permutations from the breakpoint count: `d = ⌈breakpoints / 2⌉`
(unsigned breakpoint bound, Bafna & Pevzner 1998). A breakpoint is a pair of adjacent
elements (including sentinels) whose values are not consecutive after relabelling the
target to the identity. Identical orders have 0 breakpoints and distance 0.

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L841](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L841)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `permutation1` | array of integer | Yes | Permutation 1 |
| `permutation2` | array of integer | Yes | Permutation 2 (same length) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `distance` | integer | Lower-bound number of reversals |

## Errors

| Code | Message |
|------|---------|
| 1001 | Permutation cannot be null |
| 1001 | Permutations must have the same length |

## Examples

### Example 1: Identical permutations

**User Prompt:**
> Reversal distance between [1,2,3,4] and itself.

**Expected Tool Call:**
```json
{
  "tool": "reversal_distance",
  "arguments": { "permutation1": [1, 2, 3, 4], "permutation2": [1, 2, 3, 4] }
}
```

**Response:**
```json
{ "distance": 0 }
```

### Example 2: Full reversal

**User Prompt:**
> Reversal distance between [1,2,3,4] and [4,3,2,1].

**Expected Tool Call:**
```json
{
  "tool": "reversal_distance",
  "arguments": { "permutation1": [1, 2, 3, 4], "permutation2": [4, 3, 2, 1] }
}
```

**Response:**
```json
{ "distance": 1 }
```
The reversed order has 2 breakpoints ⇒ ⌈2/2⌉ = 1 reversal.

## Performance

- **Time Complexity:** O(n).
- **Space Complexity:** O(n).

## See Also

- [detect_rearrangements](detect_rearrangements.md) — inversions/insertions/deletions
- [find_syntenic_blocks](find_syntenic_blocks.md) — collinear gene runs
