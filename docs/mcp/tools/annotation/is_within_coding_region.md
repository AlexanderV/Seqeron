# is_within_coding_region

Heuristic check whether a position lies in a coding region.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `is_within_coding_region` |
| **Method ID** | `SpliceSitePredictor.IsWithinCodingRegion` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Simple heuristic that searches up to 300 nt upstream of `position` for an `AUG` start codon; if one is
found, the position is considered coding when `(position − startIndex) % 3 == frame`. Returns `false` when
no upstream start codon exists or the position is not in the requested reading frame.

## Core Documentation Reference

- Source: [SpliceSitePredictor.cs#L980](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs#L980)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA/RNA sequence |
| `position` | integer | Yes | 0-based position to test |
| `frame` | integer | No | Reading frame 0/1/2 (default 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isCoding` | boolean | true when an upstream AUG places the position in the requested frame |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Position is out of range (ArgumentOutOfRangeException) |

## Examples

### Example 1: In-frame coding position

`AUGAAAAAAAAA`, position 6, frame 0 — AUG at index 0, `(6−0) % 3 == 0`:

**Response:**
```json
{ "isCoding": true }
```

### Example 2: No upstream start codon

```json
{ "isCoding": false }
```

## Performance

- **Time Complexity:** O(1) (bounded 300-nt upstream scan)
- **Space Complexity:** O(1)

## See Also

- [predict_gene_structure](predict_gene_structure.md) — full exon/intron structure
