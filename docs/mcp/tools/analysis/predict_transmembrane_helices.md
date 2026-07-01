# predict_transmembrane_helices

Hydropathy-based transmembrane-helix prediction.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_transmembrane_helices` |
| **Method ID** | `ProteinMotifFinder.PredictTransmembraneHelices` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Predicts **transmembrane helices** from the Kyte-Doolittle hydropathy profile: a run
of sliding windows whose mean hydropathy is ≥ `threshold` forms a TM segment. Each
segment reports the start and end residue (every residue covered by an above-threshold
window) and the peak window score.

## Core Documentation Reference

- Source: [ProteinMotifFinder.cs#L739](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs#L739)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `windowSize` | integer | No | Sliding window size (default 19, ≥ 1) |
| `threshold` | number | No | Hydropathy threshold (default 1.6) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | TM segments: `{ start, end, score }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |

## Examples

### Example 1: Single poly-Leu TM stretch

**User Prompt:**
> Predict TM helices in a D₁₀-L₂₀-D₁₀ sequence.

**Expected Tool Call:**
```json
{
  "tool": "predict_transmembrane_helices",
  "arguments": { "proteinSequence": "DDDDDDDDDDLLLLLLLLLLLLLLLLLLLLDDDDDDDDDD", "windowSize": 19, "threshold": 1.6 }
}
```

**Response:**
```json
{ "items": [ { "start": 5, "end": 34, "score": 3.8 } ] }
```
The internal poly-Leu stretch yields one TM segment spanning residues 5–34; the peak
window score is the KD value for Leu (3.8).

### Example 2: No hydrophobic segment

**User Prompt:**
> Predict TM helices in a poly-Asp peptide.

**Expected Tool Call:**
```json
{
  "tool": "predict_transmembrane_helices",
  "arguments": { "proteinSequence": "DDDDDDDDDDDDDDDDDDDD", "windowSize": 19, "threshold": 1.6 }
}
```

**Response:**
```json
{ "items": [] }
```

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(number of segments).

## See Also

- [hydrophobicity_profile](hydrophobicity_profile.md) — the underlying KD profile
- [predict_signal_peptide](predict_signal_peptide.md) — signal peptide cleavage site
