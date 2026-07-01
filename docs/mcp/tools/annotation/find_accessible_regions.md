# find_accessible_regions

Identify accessible chromatin regions (ATAC-seq-like peaks) from per-position signal.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Annotation |
| **Tool Name** | `find_accessible_regions` |
| **Method ID** | `EpigeneticsAnalyzer.FindAccessibleRegions` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Scans position-sorted accessibility signal and merges runs of positions whose signal is at or above
`threshold` into peaks, splitting whenever the gap between consecutive accessible positions exceeds
`maxGap`. A peak is reported only if its span is at least `minWidth` bp. Each region carries its maximum
signal (`accessibilityScore`) and a descriptive `peakType`: **Strong** (> 0.8), **Moderate** (> 0.5), or
**Weak**. `nearbyGenes` is populated by downstream annotation and is empty here.

Defaults: `threshold = 0.5`, `minWidth = 100`, `maxGap = 50`.

## Core Documentation Reference

- Source: [EpigeneticsAnalyzer.cs#L994](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs#L994)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `accessibilitySignal` | array | Yes | Per-position `{ position, signal }` |
| `threshold` | number | No | Signal threshold for an accessible position (default 0.5) |
| `minWidth` | integer | No | Minimum region width in bp (default 100) |
| `maxGap` | integer | No | Maximum gap before splitting a region (default 50) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `regions[].start` | integer | Region start |
| `regions[].end` | integer | Region end |
| `regions[].accessibilityScore` | number | Maximum signal in the region |
| `regions[].peakType` | string | Strong / Moderate / Weak |
| `regions[].nearbyGenes` | string[] | Nearby genes (empty here) |

## Errors

| Code | Message |
|------|---------|
| 1001 | accessibilitySignal cannot be null |

## Examples

### Example 1: Single strong peak

Positions 0..200 at signal 0.9 (gaps ≤ maxGap) form one region wider than `minWidth`:

**Response:**
```json
{ "regions": [ { "start": 0, "end": 200, "accessibilityScore": 0.9, "peakType": "Strong", "nearbyGenes": [] } ] }
```

### Example 2: All below threshold

```json
{ "regions": [] }
```

## Performance

- **Time Complexity:** O(n log n) (sort by position) for n signal points
- **Space Complexity:** O(n)

## See Also

- [predict_chromatin_state](predict_chromatin_state.md) — chromatin state from histone marks
- [annotate_histone_modifications](annotate_histone_modifications.md) — per-interval chromatin state
