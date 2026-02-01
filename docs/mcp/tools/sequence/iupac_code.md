# iupac_code

Get IUPAC ambiguity code for a set of nucleotide bases.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Sequence |
| **Tool Name** | `iupac_code` |
| **Method ID** | `IupacDnaSequence.GetIupacCode` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the IUPAC ambiguity code that represents a given set of nucleotide bases. IUPAC codes are standardized single-character representations for nucleotide ambiguity, commonly used in primer design, sequence alignment, and degenerate probe design.

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `bases` | string | Yes | Nucleotide bases to encode (e.g., "AG" for purine R) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `code` | string | IUPAC ambiguity code |
| `inputBases` | string | Input bases (normalized) |

## IUPAC Codes Reference

| Code | Bases | Mnemonic |
|------|-------|----------|
| R | A, G | puRine |
| Y | C, T | pYrimidine |
| S | G, C | Strong |
| W | A, T | Weak |
| K | G, T | Keto |
| M | A, C | aMino |
| B | C, G, T | not A |
| D | A, G, T | not C |
| H | A, C, T | not G |
| V | A, C, G | not T |
| N | A, C, G, T | aNy |

## Examples

### Example 1: Purine code

**User Prompt:**
> What IUPAC code represents A or G?

**Expected Tool Call:**
```json
{
  "tool": "iupac_code",
  "arguments": {
    "bases": "AG"
  }
}
```

**Response:**
```json
{
  "code": "R",
  "inputBases": "AG"
}
```

## See Also

- [iupac_match](iupac_match.md) - Check if codes can match
- [iupac_matches](iupac_matches.md) - Check nucleotide against code
