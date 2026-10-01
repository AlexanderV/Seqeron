# find_low_complexity_regions

Entropy-thresholded low-complexity DNA regions (per-base Shannon scan, or BBDuk's k-mer entropy masking).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_low_complexity_regions` |
| **Method ID** | `SequenceComplexity.FindLowComplexityRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Finds contiguous **low-complexity DNA regions**: a region is the union (maximal run of covered positions) of the step-1 sliding windows whose
entropy is strictly below `entropyThreshold`; overlapping flagged windows merge (BBDuk `maskLowEntropy` window-union reporting).
`method = "shannon"` (default): per-base Shannon entropy in bits (equals BBDuk with `entropyk=1 entropy=t/log₂w`).
`method = "bbduk"`: exactly the bases masked by `bbduk.sh entropy=<entropyThreshold> entropymask=t entropywindow=<windowSize> entropyk=<entropyK>`
(k-mer entropy normalised by ln(windowSize − k + 1), 0–1; BBDuk defaults window 50, k 5; BBMap 40.02, 0 mismatches on 4 500 random cases).
N / IUPAC codes are accepted; windows containing them are never flagged (BBDuk `ns() < 1`). Each region reports its bounds, length,
minimum entropy and the covered subsequence. Homopolymer and simple-repeat tracts are
the typical hits.

## Core Documentation Reference

- Source: [SequenceComplexity.cs#L602](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L602), BBDuk mode [#L560](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs#L701)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence, A/C/G/T + IUPAC codes (min length 1) |
| `windowSize` | integer | No | Window size (default 64; BBDuk's `entropywindow` default is 50) |
| `entropyThreshold` | number | No | `shannon`: bits, finite ≥ 0 (default 1.0); `bbduk`: cutoff in [0, 1] |
| `method` | string | No | `shannon` (default) or `bbduk` |
| `entropyK` | integer | No | k-mer length for `bbduk` (1–15, < windowSize; default 5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | `{ start, end, length, minEntropy, sequence }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1001 | method must be 'shannon' or 'bbduk' |
| 1001 | entropyThreshold NaN / infinite / negative (bbduk: outside [0, 1]); entropyK outside 1–15 or ≥ windowSize |

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
{ "items": [ { "start": 79, "end": 145, "length": 67, "minEntropy": 0.0 } ] }
```
The 64-nt poly-A tract (zero entropy) forms one region spanning 79–145: the union of all windows with entropy < 0.5 (first flagged window starts at 79 = "C"+19A, last at 126 = 19A+"T", covering up to 145).

### Example 1b: BBDuk mode

**Expected Tool Call:**
```json
{
  "tool": "find_low_complexity_regions",
  "arguments": { "sequence": "CGGAGCCTGTTCCTGTACCATTATCTCTTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATACCCTGAAGAGGATCTACAGATGCAAAGC", "windowSize": 50, "entropyThreshold": 0.5, "method": "bbduk", "entropyK": 5 }
}
```

**Response:**
```json
{ "items": [ { "start": 11, "end": 88, "length": 78, "minEntropy": 0.2674965262413025 } ] }
```
Identical to the bases `bbduk.sh entropy=0.5 entropymask=t` (BBMap 40.02) replaces with N.

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

- **Time Complexity:** O(n · windowSize) (`shannon`); O(n) (`bbduk`, incremental).
- **Space Complexity:** O(number of regions).

## See Also

- [windowed_complexity](windowed_complexity.md)
- [mask_low_complexity](mask_low_complexity.md)
- [find_protein_low_complexity_regions](find_protein_low_complexity_regions.md)
