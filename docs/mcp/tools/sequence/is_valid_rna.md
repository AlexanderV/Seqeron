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

## Core Documentation Reference

- Source: [SequenceExtensions.cs#L225](../../../../src/Seqeron/Seqeron.Genomics/SequenceExtensions.cs#L225)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to validate (min length: 1) |

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

## Performance

- **Time Complexity:** O(n) where n is sequence length
- **Space Complexity:** O(1)

## See Also

- [rna_validate](rna_validate.md) - Full validation with error details
- [is_valid_dna](is_valid_dna.md) - DNA validation
