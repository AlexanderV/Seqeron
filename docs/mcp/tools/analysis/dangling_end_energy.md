# dangling_end_energy

Turner 2004 dangling-end stacking energy for an RNA helix end.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dangling_end_energy` |
| **Method ID** | `RnaSecondaryStructure.GetDanglingEndEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the 5' or 3' dangling-end stacking energy (kcal/mol at 37 °C) for a single unpaired
base stacked on the end of an RNA helix, using the NNDB Turner 2004 dangling-end tables. The
lookup key is `closingBase5 + danglingBase + closingBase3`; the 3' or 5' table is chosen by the
`is3Prime` flag. Combinations absent from the table (e.g. a non-canonical closing pair) return
`0.0`.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L649](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L649)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `closingBase5` | string | Yes | Closing pair 5' base (length-1) |
| `closingBase3` | string | Yes | Closing pair 3' base (length-1) |
| `danglingBase` | string | Yes | Dangling base (length-1) |
| `is3Prime` | boolean | Yes | True for a 3' dangle, false for 5' |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Dangling-end stacking energy (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: 3' dangle on a G-C helix end

**Input:** `{ "closingBase5": "G", "closingBase3": "C", "danglingBase": "A", "is3Prime": true }`
→ key `GAC` in the 3' table → **Response:** `{ "energy": -1.1 }`

### Example 2: same context, 5' dangle

**Input:** `{ "closingBase5": "G", "closingBase3": "C", "danglingBase": "A", "is3Prime": false }`
→ key `GAC` in the 5' table → **Response:** `{ "energy": -0.5 }`

## Performance

- **Time Complexity:** O(1). **Space Complexity:** O(1).

## See Also

- [terminal_mismatch_energy](terminal_mismatch_energy.md)
- [stem_energy](stem_energy.md)
