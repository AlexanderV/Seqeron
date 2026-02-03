# melting_temperature

Calculate simple melting temperature of a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `melting_temperature` |
| **Method ID** | `SequenceStatistics.CalculateMeltingTemperature` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the melting temperature (Tm) of a DNA sequence using either the Wallace rule for short oligonucleotides or the Marmur-Doty GC formula for longer sequences. The Wallace rule uses `Tm = 2(A+T) + 4(G+C)` and is suitable for primers under 14 bp.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L441](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceStatistics.cs#L441)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence (min length: 1) |
| `useWallaceRule` | boolean | No | Use Wallace rule for short oligos (default: true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `tm` | number | Melting temperature |
| `unit` | string | Temperature unit (°C) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Short primer (Wallace rule)

**User Prompt:**
> What's the Tm of primer "ATGCGATCGATCG"?

**Expected Tool Call:**
```json
{
  "tool": "melting_temperature",
  "arguments": {
    "sequence": "ATGCGATCGATCG"
  }
}
```

**Response:**
```json
{
  "tm": 38,
  "unit": "°C"
}
```

### Example 2: Using GC formula

**User Prompt:**
> Calculate Tm using GC formula for "ATGCGATCGATCG"

**Expected Tool Call:**
```json
{
  "tool": "melting_temperature",
  "arguments": {
    "sequence": "ATGCGATCGATCG",
    "useWallaceRule": false
  }
}
```

**Response:**
```json
{
  "tm": 60.5,
  "unit": "°C"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [thermodynamics](thermodynamics.md) - Full thermodynamic analysis
- [nucleotide_composition](nucleotide_composition.md) - GC content analysis
