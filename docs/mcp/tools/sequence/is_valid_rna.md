# is_valid_rna

Quick validation if a sequence contains only valid RNA characters.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `is_valid_rna` |
| **Method ID** | `SequenceExtensions.IsValidRna` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs a fast check if a sequence contains only valid RNA characters (A, U, G, C). This is faster than `rna_validate` but returns less detailed information. Use this when you only need a boolean check without error details.

With `iupac=true` the alphabet widens to the 15 IUPAC nucleotide codes (Biopython 1.88 `IUPACData.ambiguous_rna_letters` (`GAUCRYWSMKHBVDN`) = scikit-bio `RNA` definite ∪ degenerate chars), delegating to Core `SequenceExtensions.IsValidIupacRna`; T, `X` and gap symbols remain invalid. The default (`false`) keeps the strict behaviour.

**Empty input** is rejected with error 1001, as by every Seqeron MCP tool: null/empty input → `ArgumentException` is the project-wide MCP input convention (`docs/mcp-prompt.md`, Definition of Done §1–§2). This deliberately differs from the Core predicates (`SequenceExtensions.IsValid*` return `true` for an empty span), Biopython and scikit-bio, which treat a zero-length sequence as valid.

## Core Documentation Reference

- Source: [SequenceExtensions.cs#L225](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs#L225)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to validate (min length: 1) |
| `iupac` | boolean | No | Also accept IUPAC ambiguity codes R, Y, S, W, K, M, B, D, H, V, N (T, X and gaps still invalid). Default false = strict. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isValid` | boolean | True if all characters are valid RNA |
| `length` | integer | Sequence length |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Valid RNA

**User Prompt:**
> Is "AUGCGAUCG" valid RNA?

**Expected Tool Call:**
```json
{
  "tool": "is_valid_rna",
  "arguments": {
    "sequence": "AUGCGAUCG"
  }
}
```

**Response:**
```json
{
  "isValid": true,
  "length": 9
}
```

### Example 2: Invalid RNA (contains DNA)

**User Prompt:**
> Check if "ATGC" is valid RNA

**Expected Tool Call:**
```json
{
  "tool": "is_valid_rna",
  "arguments": {
    "sequence": "ATGC"
  }
}
```

**Response:**
```json
{
  "isValid": false,
  "length": 4
}
```

### Example 3: IUPAC ambiguity codes accepted with `iupac=true`

**Expected Tool Call:**
```json
{
  "tool": "is_valid_rna",
  "arguments": {
    "sequence": "ACGURYN",
    "iupac": true
  }
}
```

**Response:**
```json
{
  "isValid": true,
  "length": 7
}
```

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [rna_validate](rna_validate.md) - Full validation with error details
- [is_valid_dna](is_valid_dna.md) - DNA validation
