# predict_chou_fasman

Per-window Chou-Fasman helix/sheet/turn propensities for a protein.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `predict_chou_fasman` |
| **Method ID** | `SequenceStatistics.PredictSecondaryStructure` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the mean **Chou-Fasman conformational propensities** — helix (Pa), sheet
(Pb), and turn (Pt) — over each sliding window of `windowSize` residues (Chou &
Fasman, 1978). A window mean above 1.0 indicates a residue stretch that favours the
corresponding conformation. Unknown residues (X, B, Z, gaps) are skipped and excluded
from the window average. When the window exceeds the sequence length the result is
empty.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L840](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L840)
- Evidence: `docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence, one-letter code (min length 1) |
| `windowSize` | integer | No | Sliding window size (default 7, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | One `(helix, sheet, turn)` mean-propensity triple per window position |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |

## Examples

### Example 1: Single-residue window (Alanine)

**User Prompt:**
> Chou-Fasman propensities of "A" with window 1.

**Expected Tool Call:**
```json
{
  "tool": "predict_chou_fasman",
  "arguments": { "proteinSequence": "A", "windowSize": 1 }
}
```

**Response:**
```json
{ "items": [ { "helix": 1.42, "sheet": 0.83, "turn": 0.66 } ] }
```
Alanine's published Chou-Fasman parameters (Pa, Pb, Pt).

### Example 2: Two windows (AE, window 1)

**User Prompt:**
> Chou-Fasman propensities of "AE" residue-by-residue.

**Expected Tool Call:**
```json
{
  "tool": "predict_chou_fasman",
  "arguments": { "proteinSequence": "AE", "windowSize": 1 }
}
```

**Response:**
```json
{ "items": [ { "helix": 1.42, "sheet": 0.83, "turn": 0.66 }, { "helix": 1.51, "sheet": 0.37, "turn": 0.74 } ] }
```
Glutamate (E) is a strong helix former (Pa = 1.51).

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(number of windows).

## See Also

- [hydrophobicity_profile](hydrophobicity_profile.md) — Kyte-Doolittle hydropathy
- [predict_disorder](predict_disorder.md) — intrinsic disorder prediction
