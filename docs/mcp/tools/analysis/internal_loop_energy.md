# internal_loop_energy

Turner 2004 internal-loop free energy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `internal_loop_energy` |
| **Method ID** | `RnaSecondaryStructure.CalculateInternalLoopEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **internal-loop free energy** (kcal/mol) under the Turner 2004 model. When
both sides have exactly one unpaired base (`n1 = n2 = 1`) the 1×1 `int11` lookup table
is used; otherwise the energy is size initiation + asymmetry penalty + terminal
mismatches at the two closing pairs.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L759](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L759)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `n1` | integer | Yes | Unpaired bases on the 5' side |
| `n2` | integer | Yes | Unpaired bases on the 3' side |
| `closingBase5_1`, `closingBase3_1` | string | Yes | Outer closing pair |
| `closingBase5_2`, `closingBase3_2` | string | Yes | Inner closing pair |
| `mismatch5_1`, `mismatch3_1` | string | Yes | Unpaired bases adjacent to the outer pair |
| `mismatch5_2`, `mismatch3_2` | string | Yes | Unpaired bases adjacent to the inner pair |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Internal-loop energy (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: 1×1 CG/CG with G-G mismatch

**User Prompt:**
> Internal-loop energy for a 1×1 CG/CG loop with a G-G mismatch.

**Expected Tool Call:**
```json
{ "tool": "internal_loop_energy", "arguments": { "n1": 1, "n2": 1, "closingBase5_1": "C", "closingBase3_1": "G", "closingBase5_2": "C", "closingBase3_2": "G", "mismatch5_1": "G", "mismatch3_1": "G", "mismatch5_2": "G", "mismatch3_2": "G" } }
```

**Response:**
```json
{ "energy": -2.2 }
```
Read from the `int11` table (strongly stabilizing).

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [hairpin_loop_energy](hairpin_loop_energy.md)
- [bulge_loop_energy](bulge_loop_energy.md)
- [multibranch_loop_energy](multibranch_loop_energy.md)
