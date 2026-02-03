# thermodynamics

Calculate thermodynamic properties of a DNA duplex.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `thermodynamics` |
| **Method ID** | `SequenceStatistics.CalculateThermodynamics` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates thermodynamic properties of a DNA duplex using the nearest-neighbor method. Returns enthalpy (ΔH), entropy (ΔS), Gibbs free energy (ΔG), and melting temperature (Tm). Includes salt correction for Na+ concentration.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L381](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceStatistics.cs#L381)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence (min length: 1) |
| `naConcentration` | number | No | Na+ concentration in M (default: 0.05 = 50mM) |
| `primerConcentration` | number | No | Primer concentration in M (default: 0.00000025 = 250nM) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `deltaH` | number | Enthalpy change (kcal/mol) |
| `deltaS` | number | Entropy change (cal/mol·K) |
| `deltaG` | number | Gibbs free energy change at 37°C (kcal/mol) |
| `meltingTemperature` | number | Melting temperature (°C) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA duplex thermodynamics

**User Prompt:**
> What are the thermodynamic properties of "ATGCGATCGATCG"?

**Expected Tool Call:**
```json
{
  "tool": "thermodynamics",
  "arguments": {
    "sequence": "ATGCGATCGATCG"
  }
}
```

**Response:**
```json
{
  "deltaH": -98.5,
  "deltaS": -267.2,
  "deltaG": -15.6,
  "meltingTemperature": 42.3
}
```

### Example 2: Custom salt concentration

**User Prompt:**
> Calculate Tm for "GCGCGCGC" at 100mM Na+

**Expected Tool Call:**
```json
{
  "tool": "thermodynamics",
  "arguments": {
    "sequence": "GCGCGCGC",
    "naConcentration": 0.1
  }
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [melting_temperature](melting_temperature.md) - Simple Tm calculation
- [nucleotide_composition](nucleotide_composition.md) - GC content analysis
