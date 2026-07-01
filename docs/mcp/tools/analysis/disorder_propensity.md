# disorder_propensity

TOP-IDP intrinsic-disorder propensity of a single amino acid.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `disorder_propensity` |
| **Method ID** | `DisorderPredictor.GetDisorderPropensity` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the TOP-IDP disorder-propensity value for a single amino acid (Campen et al. 2008,
Table 2). Higher values indicate greater propensity for intrinsic disorder; the scale runs from
`W = −0.884` (most order-promoting) to `P = 0.987` (most disorder-promoting). The lookup is
case-insensitive; an unrecognized residue returns `0`.

## Core Documentation Reference

- Source: [DisorderPredictor.cs#L865](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs#L865)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `aminoAcid` | string | Yes | Single amino-acid letter (length-1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `propensity` | number | TOP-IDP propensity value |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: Proline (most disorder-promoting)

**Input:** `{ "aminoAcid": "P" }` → **Response:** `{ "propensity": 0.987 }`

### Example 2: Tryptophan (most order-promoting)

**Input:** `{ "aminoAcid": "W" }` → **Response:** `{ "propensity": -0.884 }`

## Performance

- **Time Complexity:** O(1). **Space Complexity:** O(1).

## See Also

- [is_disorder_promoting](is_disorder_promoting.md)
- [predict_disorder](predict_disorder.md)
