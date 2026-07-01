# predict_chromatin_state

Predict chromatin state from histone-modification signal levels.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `predict_chromatin_state` |
| **Method ID** | `EpigeneticsAnalyzer.PredictChromatinState` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Applies a ChromHMM-style present/absent rule set (Ernst & Kellis 2012; Roadmap Epigenomics) to six histone
marks. Each mark is binarised at a presence threshold of 0.5 and the state is decided in priority order:

1. H3K4me3 + H3K27me3 → **BivalentPromoter**
2. H3K4me1 + H3K27me3 → **BivalentEnhancer**
3. H3K4me3 → **ActivePromoter**
4. H3K4me1 → **ActiveEnhancer** (with H3K27ac) or **WeakEnhancer** (without)
5. H3K36me3 → **Transcribed**
6. H3K27me3 → **Repressed**
7. H3K9me3 → **Heterochromatin**
8. none → **LowSignal**

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L886](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L886)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `h3k4me3` | number | Yes | H3K4me3 signal (0..1) |
| `h3k4me1` | number | Yes | H3K4me1 signal (0..1) |
| `h3k27ac` | number | Yes | H3K27ac signal (0..1) |
| `h3k36me3` | number | Yes | H3K36me3 signal (0..1) |
| `h3k27me3` | number | Yes | H3K27me3 signal (0..1) |
| `h3k9me3` | number | Yes | H3K9me3 signal (0..1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `state` | string | Predicted chromatin state (see priority list) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Histone signal must be in [0, 1] (ArgumentOutOfRangeException) |

## Examples

### Example 1: Active promoter

**User Prompt:**
> Only H3K4me3 is high — what chromatin state is this?

**Response:**
```json
{ "state": "ActivePromoter" }
```

### Example 2: Bivalent promoter

With both H3K4me3 and H3K27me3 present, the bivalent rule wins:

**Response:**
```json
{ "state": "BivalentPromoter" }
```

## Performance

- **Time Complexity:** O(1)
- **Space Complexity:** O(1)

## See Also

- [annotate_histone_modifications](annotate_histone_modifications.md) — per-interval chromatin state from a single mark
- [find_accessible_regions](find_accessible_regions.md) — accessible chromatin peaks
