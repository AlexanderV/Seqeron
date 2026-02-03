# molecular_weight_nucleotide

Calculate molecular weight of a DNA or RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `molecular_weight_nucleotide` |
| **Method ID** | `SequenceStatistics.CalculateNucleotideMolecularWeight` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the molecular weight of a DNA or RNA sequence in Daltons (Da). Uses average molecular weights for nucleotides as monophosphates. Accounts for different weights between DNA (deoxyribonucleotides) and RNA (ribonucleotides).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L180](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceStatistics.cs#L180)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA or RNA sequence (min length: 1) |
| `isDna` | boolean | No | True for DNA, false for RNA (default: true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `molecularWeight` | number | Molecular weight in Daltons |
| `unit` | string | Unit of measurement ("Da") |
| `sequenceType` | string | "DNA" or "RNA" |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: DNA molecular weight

**User Prompt:**
> What's the molecular weight of DNA "ATGC"?

**Expected Tool Call:**
```json
{
  "tool": "molecular_weight_nucleotide",
  "arguments": {
    "sequence": "ATGC",
    "isDna": true
  }
}
```

**Response:**
```json
{
  "molecularWeight": 1307.8,
  "unit": "Da",
  "sequenceType": "DNA"
}
```

### Example 2: RNA molecular weight

**User Prompt:**
> Calculate molecular weight of RNA "AUGC"

**Expected Tool Call:**
```json
{
  "tool": "molecular_weight_nucleotide",
  "arguments": {
    "sequence": "AUGC",
    "isDna": false
  }
}
```

**Response:**
```json
{
  "molecularWeight": 1357.2,
  "unit": "Da",
  "sequenceType": "RNA"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [molecular_weight_protein](molecular_weight_protein.md) - Protein molecular weight
- [nucleotide_composition](nucleotide_composition.md) - DNA/RNA composition
