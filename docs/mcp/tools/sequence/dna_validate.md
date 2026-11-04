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

With `iupac=true` the alphabet widens to the 15 IUPAC nucleotide codes (Biopython 1.88 `IUPACData.ambiguous_dna_letters` (`GATCRYWSMKHBVDN`) = scikit-bio `DNA` definite ∪ degenerate chars), delegating to Core `SequenceExtensions.IndexOfInvalidIupacDna`; U, `X` and gap symbols remain invalid. The default (`false`) keeps the strict behaviour.

**Empty input** is rejected with error 1001, as by every Seqeron MCP tool: null/empty input → `ArgumentException` is the project-wide MCP input convention (`docs/mcp-prompt.md`, Definition of Done §1–§2). This deliberately differs from the Core predicates (`SequenceExtensions.IsValid*` return `true` for an empty span), Biopython and scikit-bio, which treat a zero-length sequence as valid.

## Core Documentation Reference

- Source: [DnaSequence.cs#L115](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs#L115)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The DNA sequence to validate (min length: 1) |
| `iupac` | boolean | No | Also accept IUPAC ambiguity codes R, Y, S, W, K, M, B, D, H, V, N (U, X and gaps still invalid). Default false = strict A/C/G/T. |

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

### Example 3: IUPAC ambiguity codes accepted with `iupac=true`

**Expected Tool Call:**
```json
{
  "tool": "dna_validate",
  "arguments": {
    "sequence": "ACGTRYN",
    "iupac": true
  }
}
```

**Response:**
```json
{
  "valid": true,
  "length": 7,
  "error": null
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [rna_validate](rna_validate.md) - Validate RNA sequences
- [dna_reverse_complement](dna_reverse_complement.md) - Get reverse complement
