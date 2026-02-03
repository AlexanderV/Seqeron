# protein_validate

Validate a protein (amino acid) sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `protein_validate` |
| **Method ID** | `ProteinSequence.TryCreate` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Validates whether a sequence contains only valid amino acid codes. Valid characters include the 20 standard amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y), stop codon (*), and unknown (X). Case-insensitive validation.

## Core Documentation Reference

- Source: [ProteinSequence.cs#L357](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/ProteinSequence.cs#L357)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The protein sequence to validate (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `valid` | boolean | Whether the sequence is a valid protein |
| `length` | integer | Length of the sequence |
| `error` | string? | Error message if invalid, null if valid |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Valid protein sequence

**User Prompt:**
> Is "MAEGEITTFT" a valid protein?

**Expected Tool Call:**
```json
{
  "tool": "protein_validate",
  "arguments": {
    "sequence": "MAEGEITTFT"
  }
}
```

**Response:**
```json
{
  "valid": true,
  "length": 10,
  "error": null
}
```

### Example 2: Invalid protein sequence

**User Prompt:**
> Validate protein "MAEGEITTFJ"

**Expected Tool Call:**
```json
{
  "tool": "protein_validate",
  "arguments": {
    "sequence": "MAEGEITTFJ"
  }
}
```

**Response:**
```json
{
  "valid": false,
  "length": 10,
  "error": "Invalid amino acid 'J' at position 9"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [amino_acid_composition](amino_acid_composition.md) - Get amino acid frequencies
- [dna_validate](dna_validate.md) - Validate DNA sequences
