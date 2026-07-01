# annotate_histone_modifications

Annotate intervals with predicted chromatin state from a histone mark and signal level.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `annotate_histone_modifications` |
| **Method ID** | `EpigeneticsAnalyzer.AnnotateHistoneModifications` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Labels each input interval with the canonical Roadmap Epigenomics chromatin state implied by its
single histone mark. A mark with signal below the presence threshold (`0.5`, ChromHMM-style
binarization; Ernst & Kellis 2012) is annotated `LowSignal`. Otherwise the mark maps to its
canonical state:

| Mark | State |
|------|-------|
| H3K4me3 | ActivePromoter (TssA) |
| H3K4me1 | WeakEnhancer (Enh, without H3K27ac) |
| H3K27ac | ActiveEnhancer |
| H3K36me3 | Transcribed (Tx) |
| H3K27me3 | Repressed (ReprPC) |
| H3K9me3 | Heterochromatin (Het) |
| H3K9ac | ActivePromoter |
| (other) | LowSignal |

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L950](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L950)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `modifications` | array | Yes | List of `{ start, end, mark, signal }` intervals |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `annotations` | array | Per-interval `{ start, end, mark, signal, predictedState }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | modifications cannot be null |

## Examples

### Example 1: Active promoter mark

**Expected Tool Call:**
```json
{
  "tool": "annotate_histone_modifications",
  "arguments": { "modifications": [ { "start": 0, "end": 100, "mark": "H3K4me3", "signal": 0.9 } ] }
}
```

**Response:**
```json
{ "annotations": [ { "start": 0, "end": 100, "mark": "H3K4me3", "signal": 0.9, "predictedState": "ActivePromoter" } ] }
```

### Example 2: Below-threshold mark is LowSignal

**Response:**
```json
{ "annotations": [ { "start": 0, "end": 100, "mark": "H3K4me3", "signal": 0.3, "predictedState": "LowSignal" } ] }
```

## Performance

- **Time Complexity:** O(n) in number of intervals
- **Space Complexity:** O(n)

## See Also

- [predict_chromatin_state](predict_chromatin_state.md) - Combine multiple marks into one state
- [find_accessible_regions](find_accessible_regions.md) - ATAC-seq-like accessible peaks
