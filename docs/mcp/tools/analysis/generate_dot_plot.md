# generate_dot_plot

Matching k-mer coordinates between two sequences (dot plot).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `generate_dot_plot` |
| **Method ID** | `ComparativeGenomics.GenerateDotPlot` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Generates the coordinates for a **dot plot**: for every shared `wordSize`-mer, a point
`(i, j)` where `i` is the offset in sequence 1 and `j` an occurrence offset in sequence
2. Words are slid along sequence 1 by `stepSize` and located in sequence 2 via a suffix
tree. Diagonal runs of points reveal collinear similarity.

## Core Documentation Reference

- Source: [ComparativeGenomics.cs#L1262](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs#L1262)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | Sequence 1 (min length 1) |
| `sequence2` | string | Yes | Sequence 2 (min length 1) |
| `wordSize` | integer | No | Word size (default 10, ≥ 1) |
| `stepSize` | integer | No | Step size (default 1, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `points` | array | Match coordinates `{ x, y }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Repeated word matches

**User Prompt:**
> Dot plot of "ATGCATGC" against itself, word 4, step 4.

**Expected Tool Call:**
```json
{
  "tool": "generate_dot_plot",
  "arguments": { "sequence1": "ATGCATGC", "sequence2": "ATGCATGC", "wordSize": 4, "stepSize": 4 }
}
```

**Response:**
```json
{ "points": [ { "x": 0, "y": 0 }, { "x": 0, "y": 4 }, { "x": 4, "y": 0 }, { "x": 4, "y": 4 } ] }
```
ATGC occurs at offsets 0 and 4 in both sequences.

### Example 2: No shared word

**User Prompt:**
> Dot plot of "ATGC" against "TTTT", word 4.

**Expected Tool Call:**
```json
{
  "tool": "generate_dot_plot",
  "arguments": { "sequence1": "ATGC", "sequence2": "TTTT", "wordSize": 4, "stepSize": 1 }
}
```

**Response:**
```json
{ "points": [] }
```

## Performance

- **Time Complexity:** O(n) suffix-tree build + O((m/step)·(word + occ)).
- **Space Complexity:** O(number of matches).

## See Also

- [find_syntenic_blocks](find_syntenic_blocks.md)
- [calculate_ani](calculate_ani.md)
