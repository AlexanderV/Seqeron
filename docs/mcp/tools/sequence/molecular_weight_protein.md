# molecular_weight_protein

Calculate molecular weight of a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `molecular_weight_protein` |
| **Method ID** | `SequenceStatistics.CalculateMolecularWeight` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the molecular weight of a protein sequence in Daltons (Da). Uses average isotopic masses for amino acids and accounts for water loss during peptide bond formation.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L159](../../../../src/Seqeron/Seqeron.Genomics/SequenceStatistics.cs#L159)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `molecularWeight` | number | Molecular weight in Daltons |
| `unit` | string | Unit of measurement ("Da") |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2002 | Invalid protein sequence |

## Examples

### Example 1: Protein molecular weight

**User Prompt:**
> What's the molecular weight of "MAEGEITTFT"?

**Expected Tool Call:**
```json
{
  "tool": "molecular_weight_protein",
  "arguments": {
    "sequence": "MAEGEITTFT"
  }
}
```

**Response:**
```json
{
  "molecularWeight": 1140.25,
  "unit": "Da"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [molecular_weight_nucleotide](molecular_weight_nucleotide.md) - DNA/RNA molecular weight
- [amino_acid_composition](amino_acid_composition.md) - Full protein analysis
