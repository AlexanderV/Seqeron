# multibranch_loop_energy

Turner 2004 multibranch-loop free energy.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `multibranch_loop_energy` |
| **Method ID** | `RnaSecondaryStructure.CalculateMultibranchLoopEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **multibranch-loop free energy** (kcal/mol) with the Turner 2004 affine
model: `a + b·(unpaired/helices) + c·helices`, where `a = 9.25`, `b = 0.91`,
`c = −0.63`. An optional pre-computed stacking/dangling term is added, and a strain
penalty of +3.14 is applied when the junction is sterically strained.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L961](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L961)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `numHelices` | integer | Yes | Number of helical branches |
| `numUnpaired` | integer | Yes | Total unpaired nucleotides in the loop |
| `hasStrain` | boolean | No | Junction has steric strain (default false) |
| `stackingEnergy` | number | No | Pre-computed optimal stacking/dangling energy (default 0) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Multibranch-loop energy (kcal/mol) |

## Errors

None (numeric inputs).

## Examples

### Example 1: 3-way junction, 6 unpaired

**User Prompt:**
> Multibranch energy for a 3-way junction with 6 unpaired bases.

**Expected Tool Call:**
```json
{ "tool": "multibranch_loop_energy", "arguments": { "numHelices": 3, "numUnpaired": 6 } }
```

**Response:**
```json
{ "energy": 9.18 }
```
9.25 + 0.91·(6/3) + (−0.63)·3 = 9.18.

### Example 2: Strained junction with stacking

**User Prompt:**
> Multibranch energy for a strained 3-way junction with 1 unpaired base and −2.0 stacking.

**Expected Tool Call:**
```json
{ "tool": "multibranch_loop_energy", "arguments": { "numHelices": 3, "numUnpaired": 1, "hasStrain": true, "stackingEnergy": -2.0 } }
```

**Response:**
```json
{ "energy": 8.8 }
```
Adds the −2.0 stacking term and the +3.14 strain penalty.

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [flush_coaxial_stacking](flush_coaxial_stacking.md)
- [mismatch_coaxial_stacking](mismatch_coaxial_stacking.md)
