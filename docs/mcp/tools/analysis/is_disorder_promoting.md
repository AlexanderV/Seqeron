# is_disorder_promoting

Whether an amino acid is in Dunker's disorder-promoting set.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `is_disorder_promoting` |
| **Method ID** | `DisorderPredictor.IsDisorderPromoting` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns whether an amino acid belongs to Dunker's disorder-promoting set `{A, R, G, Q, S, P, E, K}`
(Dunker et al. 2001). Order-promoting residues `{W, C, F, I, Y, V, L, N}` and ambiguous residues
`{H, M, T, D}` return `false`. The lookup is case-insensitive; unrecognized residues return `false`.

## Core Documentation Reference

- Source: [DisorderPredictor.cs#L876](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs#L876)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `aminoAcid` | string | Yes | Single amino-acid letter (length-1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `result` | boolean | True if the residue is disorder-promoting |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: Proline (disorder-promoting)

**Input:** `{ "aminoAcid": "P" }` → **Response:** `{ "result": true }`

### Example 2: Tryptophan (order-promoting)

**Input:** `{ "aminoAcid": "W" }` → **Response:** `{ "result": false }`

## Performance

- **Time Complexity:** O(1). **Space Complexity:** O(1).

## See Also

- [disorder_propensity](disorder_propensity.md)
- [predict_disorder](predict_disorder.md)
