# isoelectric_point

Calculate isoelectric point (pI) of a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `isoelectric_point` |
| **Method ID** | `SequenceStatistics.CalculateIsoelectricPoint` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the isoelectric point (pI) of a protein sequence. The pI is the pH at which the protein has no net electrical charge. Uses binary search to find the pH where net charge equals zero, considering ionizable side chains (D, E, C, Y, H, K, R) and terminal amino/carboxyl groups.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L228](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceStatistics.cs#L228)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pI` | number | Isoelectric point (pH value 0-14) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2002 | Invalid protein sequence |

## Examples

### Example 1: Acidic protein

**User Prompt:**
> What's the isoelectric point of "MAEGEITTFT"?

**Expected Tool Call:**
```json
{
  "tool": "isoelectric_point",
  "arguments": {
    "sequence": "MAEGEITTFT"
  }
}
```

**Response:**
```json
{
  "pI": 4.21
}
```

### Example 2: Basic protein

**User Prompt:**
> Calculate pI of "MKRLLKR"

**Expected Tool Call:**
```json
{
  "tool": "isoelectric_point",
  "arguments": {
    "sequence": "MKRLLKR"
  }
}
```

**Response:**
```json
{
  "pI": 11.5
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [amino_acid_composition](amino_acid_composition.md) - Full protein composition analysis
- [protein_validate](protein_validate.md) - Validate protein sequences
