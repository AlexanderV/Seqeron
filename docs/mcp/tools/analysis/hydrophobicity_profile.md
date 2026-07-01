# hydrophobicity_profile

Sliding-window Kyte-Doolittle hydropathy profile of a protein.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `hydrophobicity_profile` |
| **Method ID** | `SequenceStatistics.CalculateHydrophobicityProfile` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **Kyte-Doolittle hydropathy profile**: the unweighted mean hydropathy
value over each sliding window of `windowSize` residues (Kyte & Doolittle, 1982).
Positive values indicate hydrophobic stretches (candidate membrane-spanning or buried
regions), negative values hydrophilic ones. The profile has exactly
`N − windowSize + 1` values; non-standard residues contribute 0 to a window's sum.
When the window exceeds the sequence length the result is empty.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L404](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L404)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence (min length 1) |
| `windowSize` | integer | No | Sliding window size (default 9, ≥ 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `values` | array of number | Mean Kyte-Doolittle hydropathy per window (length `N − windowSize + 1`) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Window size must be at least 1 |

## Examples

### Example 1: Hydrophobic poly-isoleucine (window 3)

**User Prompt:**
> Hydrophobicity profile of "IIIII" with window 3.

**Expected Tool Call:**
```json
{
  "tool": "hydrophobicity_profile",
  "arguments": { "proteinSequence": "IIIII", "windowSize": 3 }
}
```

**Response:**
```json
{ "values": [4.5, 4.5, 4.5] }
```
Isoleucine has KD value 4.5, so every window averages 4.5.

### Example 2: Mixed residues (window 2)

**User Prompt:**
> Hydrophobicity of "AV" with window 2.

**Expected Tool Call:**
```json
{
  "tool": "hydrophobicity_profile",
  "arguments": { "proteinSequence": "AV", "windowSize": 2 }
}
```

**Response:**
```json
{ "values": [3.0] }
```
(A = 1.8, V = 4.2) ⇒ mean = 3.0.

## Performance

- **Time Complexity:** O(n · windowSize).
- **Space Complexity:** O(n) for the profile.

## See Also

- [predict_transmembrane_helices](predict_transmembrane_helices.md) — hydropathy-based TM prediction
- [predict_chou_fasman](predict_chou_fasman.md) — secondary-structure propensities
