# amino_acid_composition

Calculate amino acid composition of a protein sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `amino_acid_composition` |
| **Method ID** | `SequenceStatistics.CalculateAminoAcidComposition` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Calculates the amino acid composition of a protein sequence, including counts of each amino acid, molecular weight, isoelectric point (pI), hydrophobicity, charged residue ratio, and aromatic residue ratio.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L98](../../../../SuffixTree.Genomics/SequenceStatistics.cs#L98)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence to analyze (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `length` | integer | Number of amino acids |
| `counts` | object | Count of each amino acid |
| `molecularWeight` | number | Molecular weight in Daltons |
| `isoelectricPoint` | number | Isoelectric point (pI), 0-14 |
| `hydrophobicity` | number | Average hydrophobicity |
| `chargedResidueRatio` | number | Ratio of charged residues (D, E, K, R, H) |
| `aromaticResidueRatio` | number | Ratio of aromatic residues (F, Y, W) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 2002 | Invalid protein sequence |

## Examples

### Example 1: Protein composition

**User Prompt:**
> What's the amino acid composition of "MAEGEITTFT"?

**Expected Tool Call:**
```json
{
  "tool": "amino_acid_composition",
  "arguments": {
    "sequence": "MAEGEITTFT"
  }
}
```

**Response:**
```json
{
  "length": 10,
  "counts": {"M": 1, "A": 1, "E": 2, "G": 1, "I": 1, "T": 3, "F": 1},
  "molecularWeight": 1140.25,
  "isoelectricPoint": 4.21,
  "hydrophobicity": 0.12,
  "chargedResidueRatio": 0.2,
  "aromaticResidueRatio": 0.1
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1) for standard amino acids

## See Also

- [protein_validate](protein_validate.md) - Validate protein sequences
- [nucleotide_composition](nucleotide_composition.md) - DNA/RNA composition
