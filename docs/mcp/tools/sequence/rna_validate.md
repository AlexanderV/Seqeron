# rna_validate

Validate an RNA sequence.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `rna_validate` |
| **Method ID** | `RnaSequence.TryCreate` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Validates whether a sequence contains only valid RNA nucleotides (A, C, G, U). Returns validation status, sequence length, and detailed error message if invalid. Case-insensitive validation. Note: T (thymine) is invalid in RNA; use U (uracil) instead.

With `iupac=true` the alphabet widens to the 15 IUPAC nucleotide codes (Biopython 1.88 `IUPACData.ambiguous_rna_letters` (`GAUCRYWSMKHBVDN`) = scikit-bio `RNA` definite ∪ degenerate chars), delegating to Core `SequenceExtensions.IndexOfInvalidIupacRna`; T, `X` and gap symbols remain invalid. The default (`false`) keeps the strict behaviour.

**Empty input** is rejected with error 1001, as by every Seqeron MCP tool: null/empty input → `ArgumentException` is the project-wide MCP input convention (`docs/mcp-prompt.md`, Definition of Done §1–§2). This deliberately differs from the Core predicates (`SequenceExtensions.IsValid*` return `true` for an empty span), Biopython and scikit-bio, which treat a zero-length sequence as valid.

## Core Documentation Reference

- Source: [RnaSequence.cs#L169](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/RnaSequence.cs#L169)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The RNA sequence to validate (min length: 1) |
| `iupac` | boolean | No | Also accept IUPAC ambiguity codes R, Y, S, W, K, M, B, D, H, V, N (T, X and gaps still invalid). Default false = strict A/C/G/U. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `valid` | boolean | Whether the sequence is valid RNA |
| `length` | integer | Length of the sequence |
| `error` | string? | Error message if invalid, null if valid |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Valid RNA sequence

**User Prompt:**
> Is "AUGCAUGC" a valid RNA sequence?

**Expected Tool Call:**
```json
{
  "tool": "rna_validate",
  "arguments": {
    "sequence": "AUGCAUGC"
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

### Example 2: Invalid RNA sequence (contains T)

**User Prompt:**
> Validate the RNA sequence "AUGTATGC"

**Expected Tool Call:**
```json
{
  "tool": "rna_validate",
  "arguments": {
    "sequence": "AUGTATGC"
  }
}
```

**Response:**
```json
{
  "valid": false,
  "length": 8,
  "error": "Invalid nucleotide 'T' at position 3"
}
```

### Example 3: IUPAC ambiguity codes accepted with `iupac=true`

**Expected Tool Call:**
```json
{
  "tool": "rna_validate",
  "arguments": {
    "sequence": "ACGURYN",
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

- [dna_validate](dna_validate.md) - Validate DNA sequences
- [rna_from_dna](rna_from_dna.md) - Transcribe DNA to RNA
