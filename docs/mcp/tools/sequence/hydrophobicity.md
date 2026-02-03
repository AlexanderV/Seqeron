# hydrophobicity

Calculate hydrophobicity (GRAVY index) of a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `hydrophobicity` |
| **Method ID** | `SequenceStatistics.CalculateHydrophobicity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the grand average of hydropathy (GRAVY) index for a protein sequence using the Kyte-Doolittle hydrophobicity scale. Positive values indicate hydrophobic proteins, while negative values indicate hydrophilic proteins.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L306](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceStatistics.cs#L306)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `gravy` | number | GRAVY index (positive=hydrophobic, negative=hydrophilic) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2002 | Invalid protein sequence |

## Examples

### Example 1: Mixed protein

**User Prompt:**
> What's the hydrophobicity of "MAEGEITTFT"?

**Expected Tool Call:**
```json
{
  "tool": "hydrophobicity",
  "arguments": {
    "sequence": "MAEGEITTFT"
  }
}
```

**Response:**
```json
{
  "gravy": -0.12
}
```

### Example 2: Hydrophobic protein

**User Prompt:**
> Calculate GRAVY for "ILVILVILVV"

**Expected Tool Call:**
```json
{
  "tool": "hydrophobicity",
  "arguments": {
    "sequence": "ILVILVILVV"
  }
}
```

**Response:**
```json
{
  "gravy": 4.2
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [amino_acid_composition](amino_acid_composition.md) - Full protein composition analysis
- [isoelectric_point](isoelectric_point.md) - Calculate protein pI
