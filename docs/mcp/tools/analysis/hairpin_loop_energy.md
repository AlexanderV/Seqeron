# hairpin_loop_energy

Turner 2004 hairpin-loop free energy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `hairpin_loop_energy` |
| **Method ID** | `RnaSecondaryStructure.CalculateHairpinLoopEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **hairpin-loop free energy** (kcal/mol) from the Turner 2004 model: loop
initiation by size plus a terminal mismatch, with special-case tri/tetra/hexaloop
bonuses, an all-C penalty, and a special G-U closure adjustment. A 3-nt loop uses
initiation only (no terminal mismatch).

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L664](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L664)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `loopSequence` | string | Yes | Loop residues between the closing pair |
| `closingBase5` | string | Yes | Closing pair 5' base (length-1) |
| `closingBase3` | string | Yes | Closing pair 3' base (length-1) |
| `specialGUClosure` | boolean | No | Apply special G-U closure bonus (default false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Hairpin-loop energy (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: 6-nt loop, A-U closure

**User Prompt:**
> Hairpin loop energy for AAAAAA closed by A-U.

**Expected Tool Call:**
```json
{ "tool": "hairpin_loop_energy", "arguments": { "loopSequence": "AAAAAA", "closingBase5": "A", "closingBase3": "U" } }
```

**Response:**
```json
{ "energy": 4.6 }
```
Initiation (5.4) + terminal mismatch AAAU (−0.8) = 4.6.

### Example 2: Special triloop

**Response:**
```json
{ "energy": 6.8 }
```

## Performance

- **Time Complexity:** O(loop length).
- **Space Complexity:** O(1).

## See Also

- [terminal_mismatch_energy](terminal_mismatch_energy.md)
- [internal_loop_energy](internal_loop_energy.md)
- [minimum_free_energy](minimum_free_energy.md)
