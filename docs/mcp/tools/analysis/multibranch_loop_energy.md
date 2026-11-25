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

Computes the **multibranch-loop free energy** (kcal/mol) in one of two sourced Turner 2004 models:

- **Linear model (default; the MFE / partition-function engines' model):**
  `9.3 + 0·numUnpaired − 0.9·numHelices + stackingEnergy` (Mathews et al. 2004; ViennaRNA 2.7.2
  `ML_closing37 = 930`, `ML_BASE37 = 0`, `ML_intern37 = −90`). `hasStrain = true` is rejected here.
- **NNDB efn2 model (when `unpairedBetweenHelices` is given):**
  `9.25 + 0.91·min(2, average asymmetry) − 0.63·numHelices + 3.14·[3-way loop with < 2 unpaired] + stackingEnergy`,
  where the average asymmetry is the mean over the helices of |unpaired 5' of the helix − unpaired 3' of the helix|
  (NNDB "Multibranch Loops" §9.1 and parameter table). `unpairedBetweenHelices` must have `numHelices` entries
  summing to `numUnpaired`, and `hasStrain` must equal the derived strain condition.

`numHelices` counts the closing helix (≥ 3).

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs) (`CalculateMultibranchLoopEnergy`, both overloads)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `numHelices` | integer | Yes | Helices in the loop, the closing helix included (≥ 3) |
| `numUnpaired` | integer | Yes | Total unpaired nucleotides in the loop (≥ 0) |
| `hasStrain` | boolean | No | efn2 strain flag; only with `unpairedBetweenHelices` (default false) |
| `stackingEnergy` | number | No | Pre-computed optimal stacking/dangling energy (default 0) |
| `unpairedBetweenHelices` | integer[] | No | efn2 input: unpaired nucleotides between consecutive helices, 5'→3' around the loop starting after the closing pair |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Multibranch-loop energy (kcal/mol) |

## Errors

`numHelices < 3`, `numUnpaired < 0`, `hasStrain` without `unpairedBetweenHelices`, or an inconsistent
`unpairedBetweenHelices` → argument error.

## Examples

### Example 1: 3-way junction, linear model

```json
{ "tool": "multibranch_loop_energy", "arguments": { "numHelices": 3, "numUnpaired": 6 } }
```

**Response:** `{ "energy": 6.6 }` — 9.3 − 0.9·3 (= ViennaRNA `eval_structure_verbose` "Multi loop" 660, dangles 0).

### Example 2: NNDB efn2 worked example (§9.3)

```json
{ "tool": "multibranch_loop_energy", "arguments": { "numHelices": 4, "numUnpaired": 8, "unpairedBetweenHelices": [2, 0, 1, 5] } }
```

**Response:** `{ "energy": 8.55 }` — average asymmetry 2.5 capped at 2.0: 9.25 + 0.91·2 − 0.63·4 (NNDB prints 8.6).

### Example 3: strained 3-way junction (efn2)

```json
{ "tool": "multibranch_loop_energy", "arguments": { "numHelices": 3, "numUnpaired": 1, "hasStrain": true, "stackingEnergy": -2.0, "unpairedBetweenHelices": [1, 0, 0] } }
```

**Response:** `{ "energy": 9.11 }` — 9.25 + 0.91·(2/3) − 1.89 − 2.0 + 3.14.

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [flush_coaxial_stacking](flush_coaxial_stacking.md)
- [mismatch_coaxial_stacking](mismatch_coaxial_stacking.md)
