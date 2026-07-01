# terminal_mismatch_energy

Turner 2004 terminal-mismatch stacking energy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `terminal_mismatch_energy` |
| **Method ID** | `RnaSecondaryStructure.GetTerminalMismatchEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **terminal-mismatch stacking energy** (kcal/mol) for a closing base pair
and the first mismatch pair stacked at a helix end, from the Turner 2004 (NNDB)
tables. Used in hairpin, internal-loop and coaxial-stacking calculations.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L639](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L639)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `closingBase5` | string | Yes | Closing pair 5' base (length-1) |
| `closingBase3` | string | Yes | Closing pair 3' base (length-1) |
| `mismatch5` | string | Yes | 5' mismatch base (length-1) |
| `mismatch3` | string | Yes | 3' mismatch base (length-1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Terminal-mismatch energy (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: GC closing, AA mismatch

**User Prompt:**
> Terminal mismatch energy for a G-C pair with an A-A mismatch.

**Expected Tool Call:**
```json
{ "tool": "terminal_mismatch_energy", "arguments": { "closingBase5": "G", "closingBase3": "C", "mismatch5": "A", "mismatch3": "A" } }
```

**Response:**
```json
{ "energy": -1.1 }
```

### Example 2: CG closing, AA mismatch

**Response:**
```json
{ "energy": -1.5 }
```

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [hairpin_loop_energy](hairpin_loop_energy.md)
- [mismatch_coaxial_stacking](mismatch_coaxial_stacking.md)
