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
helix end is `5'-closingBase5-3' / 3'-closingBase3-5'`: a 3' dangle sits 3' of `closingBase5`, a
5' dangle 5' of `closingBase3` (= ViennaRNA 2.7.2 `dangle3_37` / `dangle5_37`). Bases are
case-insensitive and T is read as U; an unknown dangling base (N, IUPAC, gap) takes the maximum
over A/C/G/U (ViennaRNA's 'N' column); a non-canonical closing pair returns `0.0`.

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
