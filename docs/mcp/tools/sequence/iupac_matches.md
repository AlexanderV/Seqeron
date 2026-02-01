# iupac_matches

Check if a nucleotide matches an IUPAC code.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `iupac_matches` |
| **Method ID** | `IupacHelper.MatchesIupac` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Determines if a specific nucleotide (A, C, G, T) matches an IUPAC ambiguity code. This is useful for validating sequences against degenerate patterns.

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `nucleotide` | string | Yes | The nucleotide to check (A, C, G, T) |
| `iupacCode` | string | Yes | The IUPAC code to match against |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | boolean | True if nucleotide matches the IUPAC code |
| `nucleotide` | string | Input nucleotide |
| `iupacCode` | string | IUPAC code checked |

## Examples

### Example 1: A matches R (purine)

```json
{
  "tool": "iupac_matches",
  "arguments": { "nucleotide": "A", "iupacCode": "R" }
}
```

**Response:** `{ "matches": true, "nucleotide": "A", "iupacCode": "R" }`

### Example 2: Any nucleotide matches N

```json
{
  "tool": "iupac_matches",
  "arguments": { "nucleotide": "A", "iupacCode": "N" }
}
```

**Response:** `{ "matches": true, "nucleotide": "A", "iupacCode": "N" }`

## See Also

- [iupac_code](iupac_code.md) - Get IUPAC code for bases
- [iupac_match](iupac_match.md) - Check if two codes can match
