# dna_validate

Validate a DNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `dna_validate` |
| **Method ID** | `DnaSequence.TryCreate` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Validates whether a sequence contains only valid DNA nucleotides (A, C, G, T). Returns validation status, sequence length, and detailed error message if invalid. Case-insensitive validation.

## Core Documentation Reference

- Source: [DnaSequence.cs#L129](../../../../src/Seqeron/Seqeron.Genomics/DnaSequence.cs#L129)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to validate (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `valid` | boolean | Whether the sequence is valid DNA |
| `length` | integer | Length of the sequence |
| `error` | string? | Error message if invalid, null if valid |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Valid DNA sequence

**User Prompt:**
> Is "ATGCATGC" a valid DNA sequence?

**Expected Tool Call:**
```json
{
  "tool": "dna_validate",
  "arguments": {
    "sequence": "ATGCATGC"
  }
}
```

**Response:**
```json
{
  "valid": true,
  "length": 8,
  "error": null
}
```

### Example 2: Invalid DNA sequence

**User Prompt:**
> Validate the sequence "ATGXATGC"

**Expected Tool Call:**
```json
{
  "tool": "dna_validate",
  "arguments": {
    "sequence": "ATGXATGC"
  }
}
```

**Response:**
```json
{
  "valid": false,
  "length": 8,
  "error": "Invalid nucleotide 'X' at position 3"
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [rna_validate](rna_validate.md) - Validate RNA sequences
- [dna_reverse_complement](dna_reverse_complement.md) - Get reverse complement
