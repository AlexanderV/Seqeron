# stem_energy

Free energy of an RNA stem.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `stem_energy` |
| **Method ID** | `RnaSecondaryStructure.CalculateStemEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the **free energy** (kcal/mol) of an RNA stem: Turner 2004 nearest-neighbour
stacking summed over consecutive base pairs, plus AU/GU terminal penalties. A stem of a
single base pair has no stacking interaction and returns 0.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L570](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L570)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | RNA sequence containing the stem |
| `basePairs` | array | Yes | Base pairs (5'→3' order): `{ position1, position2, base1, base2, type }` |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Stem free energy (kcal/mol) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Base pairs cannot be null |

## Examples

### Example 1: Single pair (no stacking)

**User Prompt:**
> Stem energy for a single G-C pair.

**Expected Tool Call:**
```json
{ "tool": "stem_energy", "arguments": { "sequence": "GAAAC", "basePairs": [ { "position1": 0, "position2": 4, "base1": "G", "base2": "C", "type": "WatsonCrick" } ] } }
```

**Response:**
```json
{ "energy": 0.0 }
```

### Example 2: GG/CC stack

**Response:**
```json
{ "energy": -3.26 }
```
Two consecutive G-C pairs stack as GG/CC (−3.26, no terminal penalty).

## Performance

- **Time Complexity:** O(number of base pairs).
- **Space Complexity:** O(1).

## See Also

- [find_stem_loops](find_stem_loops.md)
- [minimum_free_energy](minimum_free_energy.md)
