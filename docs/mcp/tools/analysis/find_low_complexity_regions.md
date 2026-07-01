# find_low_complexity_regions

Entropy-thresholded low-complexity DNA regions.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_low_complexity_regions` |
| **Method ID** | `SequenceComplexity.FindLowComplexityRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds contiguous **low-complexity DNA regions** by merging sliding windows whose
Shannon entropy falls below `entropyThreshold`. Each region reports its bounds, length,
minimum entropy and the covered subsequence. Homopolymer and simple-repeat tracts are
the typical hits.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L255](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L255)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (min length 1) |
| `windowSize` | integer | No | Window size (default 64) |
| `entropyThreshold` | number | No | Entropy threshold (default 1.0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ start, end, length, minEntropy, sequence }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |

## Examples

### Example 1: Internal poly-A region

**User Prompt:**
> Find low-complexity regions in a sequence with a central poly-A tract (window 20, threshold 0.5).

**Expected Tool Call:**
```json
{
  "tool": "find_low_complexity_regions",
  "arguments": { "sequence": "ATGCATGC…(20×) + AAAA…(64×) + ATGCATGC…(20×)", "windowSize": 20, "entropyThreshold": 0.5 }
}
```

**Response:**
```json
{ "items": [ { "start": 79, "end": 146, "minEntropy": 0.0 } ] }
```
The 64-nt poly-A tract (zero entropy) forms one region spanning 79–146.

### Example 2: High complexity

**User Prompt:**
> Low-complexity regions in a pure ATGC repeat?

**Expected Tool Call:**
```json
{
  "tool": "find_low_complexity_regions",
  "arguments": { "sequence": "ATGCATGC…(20×)", "windowSize": 20, "entropyThreshold": 0.5 }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(number of regions).

## See Also

- [windowed_complexity](windowed_complexity.md)
- [mask_low_complexity](mask_low_complexity.md)
- [find_protein_low_complexity_regions](find_protein_low_complexity_regions.md)
