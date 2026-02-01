# iupac_match

Check if two IUPAC codes can represent the same base.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `iupac_match` |
| **Method ID** | `IupacDnaSequence.CodesMatch` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Checks if two IUPAC ambiguity codes can represent the same nucleotide base. This is useful for pattern matching with degenerate sequences.

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `code1` | string | Yes | First IUPAC code (single character) |
| `code2` | string | Yes | Second IUPAC code (single character) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | boolean | True if codes can represent the same base |
| `code1` | string | First code |
| `code2` | string | Second code |

## Examples

### Example 1: R matches A

R (purine = A or G) can represent A.

```json
{
  "tool": "iupac_match",
  "arguments": { "code1": "R", "code2": "A" }
}
```

**Response:** `{ "matches": true, "code1": "R", "code2": "A" }`

### Example 2: R does not match Y

R (purine) and Y (pyrimidine) have no overlap.

```json
{
  "tool": "iupac_match",
  "arguments": { "code1": "R", "code2": "Y" }
}
```

**Response:** `{ "matches": false, "code1": "R", "code2": "Y" }`

## See Also

- [iupac_code](iupac_code.md) - Get IUPAC code for bases
- [iupac_matches](iupac_matches.md) - Check nucleotide against code
