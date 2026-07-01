# find_conserved_elements

Identify conserved genomic elements from per-position conservation scores (runs above threshold).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_conserved_elements` |
| **Method ID** | `VariantAnnotator.FindConservedElements` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Groups per-position conservation scores by chromosome, sorts by position, and reports maximal runs of
positions whose PhastCons is at least `threshold`. A run is broken when a gap of more than 10 bp
appears. Each element that spans at least `minLength` positions (`end − start + 1`) is returned with
its chromosome, inclusive `start`/`end`, and the mean PhastCons across the run as `score`.

## Core Documentation Reference

- Source: [VariantAnnotator.cs#L1222](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs#L1222)

## Input Schema

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `scores` | array | Yes | — | Per-position conservation scores (PhastCons used for thresholding) |
| `threshold` | number | No | 0.8 | PhastCons threshold for a conserved position |
| `minLength` | integer | No | 20 | Minimum element length (positions) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `elements` | array | `{ chromosome, start, end, score }` per conserved element |

## Errors

| Code | Message |
|------|---------|
| 1001 | Scores cannot be null or empty |

## Examples

### Example 1: Contiguous high-conservation run

30 contiguous positions `chr1:100..129` with PhastCons 0.9 (threshold 0.8, minLength 20) → one
element spanning [100, 129] with mean score 0.9.

**Response:**
```json
{ "elements": [ { "chromosome": "chr1", "start": 100, "end": 129, "score": 0.9 } ] }
```

### Example 2: Below threshold

Positions with PhastCons 0.5 → no conserved elements.

## Performance

- **Time Complexity:** O(n log n) (per-chromosome sort)
- **Space Complexity:** O(n)

## See Also

- [calculate_conservation](calculate_conservation.md) - Compute per-position conservation scores
- [predict_pathogenicity](predict_pathogenicity.md) - Uses conservation as evidence
