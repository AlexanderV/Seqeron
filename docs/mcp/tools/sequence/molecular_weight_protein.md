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

Calculates the molecular weight of a protein sequence in Daltons (Da). Uses average isotopic masses for amino acids and accounts for water loss during peptide bond formation. Accepts the alphabet of [protein_validate](protein_validate.md); selenocysteine (U, 168.0532 Da) and pyrrolysine (O, 255.3134 Da) are weighed with the Biopython `IUPACData.protein_weights` masses (e.g. `MKUAO` → 735.7975 Da, as Biopython `molecular_weight("MKUAO", "protein")`), while B/Z/J/X/`*` carry no mass and are skipped.

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L412](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L412)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence (min length: 1) |
| `monoisotopic` | boolean | No | Monoisotopic residue masses, water 18.010565 (default false). GAGA → 274.127717 |
| `circular` | boolean | No | Head-to-tail cyclic peptide, one more water lost (default false). GAGA → 256.2584 |
| `strict` | boolean | No | Reject B/Z/J/X/'*' (no mass) instead of skipping them (default false) |

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
| 1001 | strict: '<c>' is not a valid unambiguous letter for protein (strict mode). |

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
