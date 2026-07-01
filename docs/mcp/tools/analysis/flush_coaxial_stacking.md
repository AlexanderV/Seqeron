# flush_coaxial_stacking

Coaxial stacking energy for two flush RNA helices.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `flush_coaxial_stacking` |
| **Method ID** | `RnaSecondaryStructure.CalculateFlushCoaxialStacking` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **flush coaxial stacking energy** (kcal/mol) for two RNA helices stacked
end-to-end with no intervening unpaired bases, using the Turner 2004 Watson-Crick / G-U
stacking table.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L980](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L980)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `base5_1` | string | Yes | 5' base of first helix end |
| `base3_1` | string | Yes | 3' base of first helix end |
| `base5_2` | string | Yes | 5' base of second helix end |
| `base3_2` | string | Yes | 3' base of second helix end |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Coaxial stacking energy (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: GC onto CG

**User Prompt:**
> Flush coaxial stacking energy of GC onto CG.

**Expected Tool Call:**
```json
{ "tool": "flush_coaxial_stacking", "arguments": { "base5_1": "G", "base3_1": "C", "base5_2": "C", "base3_2": "G" } }
```

**Response:**
```json
{ "energy": -3.42 }
```
The most stable Watson-Crick stack (GC/CG).

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [mismatch_coaxial_stacking](mismatch_coaxial_stacking.md)
- [multibranch_loop_energy](multibranch_loop_energy.md)
