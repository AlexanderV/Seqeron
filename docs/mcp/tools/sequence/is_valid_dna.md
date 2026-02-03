# is_valid_dna

Quick validation if a sequence contains only valid DNA characters.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `is_valid_dna` |
| **Method ID** | `SequenceExtensions.IsValidDna` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Performs a fast check if a sequence contains only valid DNA characters (A, T, G, C). This is faster than `dna_validate` but returns less detailed information. Use this when you only need a boolean check without error details.

## Core Documentation Reference

- Source: [SequenceExtensions.cs#L210](../../../../src/Seqeron/Algorithms/Seqeron.Genomics/SequenceExtensions.cs#L210)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | The sequence to validate (min length: 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isValid` | boolean | True if all characters are valid DNA |
| `length` | integer | Sequence length |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |

## Examples

### Example 1: Valid DNA

**User Prompt:**
> Is "ATGCGATCG" valid DNA?

**Expected Tool Call:**
```json
{
  "tool": "is_valid_dna",
  "arguments": {
    "sequence": "ATGCGATCG"
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

### Example 2: Invalid DNA (contains RNA)

**User Prompt:**
> Check if "AUGC" is valid DNA

**Expected Tool Call:**
```json
{
  "tool": "is_valid_dna",
  "arguments": {
    "sequence": "AUGC"
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

- [dna_validate](dna_validate.md) - Full validation with error details
- [is_valid_rna](is_valid_rna.md) - RNA validation
