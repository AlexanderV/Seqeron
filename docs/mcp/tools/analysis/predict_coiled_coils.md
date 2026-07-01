# predict_coiled_coils

Heptad-repeat-based coiled-coil prediction.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_coiled_coils` |
| **Method ID** | `ProteinMotifFinder.PredictCoiledCoils` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts **coiled-coil regions** from the heptad (abcdefg) repeat periodicity: sliding
windows whose heptad score is ≥ `threshold` are reported as coiled-coil segments with
their start, end and peak score. A window shorter than `windowSize` cannot be scored,
so sequences below the window length yield no regions.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L976](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L976)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `windowSize` | integer | No | Sliding window size (default 28, ≥ 1) |
| `threshold` | number | No | Score threshold (default 0.5) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | Coiled-coil regions: `{ start, end, score }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |

## Examples

### Example 1: Perfect heptad repeat

**User Prompt:**
> Predict coiled coils in (LAALAAA)₅.

**Expected Tool Call:**
```json
{
  "tool": "predict_coiled_coils",
  "arguments": { "proteinSequence": "LAALAAALAALAAALAALAAALAALAAALAALAAA", "windowSize": 28, "threshold": 0.5 }
}
```

**Response:**
```json
{ "items": [ { "start": 0, "end": 34, "score": 1.0 } ] }
```
A perfect heptad periodicity scores 1.0 over the whole 35-residue span.

### Example 2: No coiled coil

**User Prompt:**
> Predict coiled coils in a poly-glycine peptide.

**Expected Tool Call:**
```json
{
  "tool": "predict_coiled_coils",
  "arguments": { "proteinSequence": "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "windowSize": 28, "threshold": 0.5 }
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

- [predict_transmembrane_helices](predict_transmembrane_helices.md) — TM helix prediction
- [find_protein_motifs](find_protein_motifs.md) — leucine-zipper and other motifs
