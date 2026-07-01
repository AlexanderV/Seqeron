# mismatch_coaxial_stacking

Mismatch-mediated coaxial stacking energy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mismatch_coaxial_stacking` |
| **Method ID** | `RnaSecondaryStructure.CalculateMismatchCoaxialStacking` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Returns the **mismatch-mediated coaxial stacking energy** (kcal/mol): a terminal
mismatch term plus a base contribution of −2.1 plus a closing-pair bonus (−0.4 for
Watson-Crick, −0.2 for G-U), per the Turner 2004 coaxial-stacking model.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L992](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L992)

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
| `energy` | number | Coaxial stacking energy (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: GC closing, AA mismatch

**User Prompt:**
> Mismatch coaxial stacking for a G-C closing pair with an A-A mismatch.

**Expected Tool Call:**
```json
{ "tool": "mismatch_coaxial_stacking", "arguments": { "closingBase5": "G", "closingBase3": "C", "mismatch5": "A", "mismatch3": "A" } }
```

**Response:**
```json
{ "energy": -3.6 }
```
tm(GAAC) −1.1 + base −2.1 + WC bonus −0.4 = −3.6.

### Example 2: GU closing, AA mismatch

**Response:**
```json
{ "energy": -2.6 }
```
tm(GAAU) −0.3 + base −2.1 + GU bonus −0.2 = −2.6.

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [terminal_mismatch_energy](terminal_mismatch_energy.md)
- [flush_coaxial_stacking](flush_coaxial_stacking.md)
