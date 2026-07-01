# bulge_loop_energy

Free energy of an RNA bulge loop (Turner 2004).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `bulge_loop_energy` |
| **Method ID** | `RnaSecondaryStructure.CalculateBulgeLoopEnergy` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Computes the free energy (kcal/mol at 37 °C) of an RNA bulge loop using the Turner
2004 nearest-neighbour parameters (NNDB `turner04/bulge.html`):

- **n = 1:** `initiation(1) + stacking(as if no bulge) + special-C bonus − RT·ln(numStates)`.
  The special-C bonus (−0.9) applies when the bulged base is C adjacent to at least
  one paired C; the degeneracy term applies only when `numStates > 1`.
- **n > 1:** `initiation(n) + terminal AU/GU penalty on each closing pair that is A-U or G-U`.

The result is rounded to two decimals.

## Core Documentation Reference

- Source: [RnaSecondaryStructure.cs#L894](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs#L894)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `bulgeSize` | integer | Yes | Number of unpaired bases in the bulge |
| `bulgedBase` | string | Yes | The bulged nucleotide (single char; used for special-C check when n=1) |
| `pair5_base1` | string | Yes | 5' base of the pair on the 5' side of the bulge |
| `pair5_base2` | string | Yes | 3' base of the pair on the 5' side |
| `pair3_base1` | string | Yes | 5' base of the pair on the 3' side of the bulge |
| `pair3_base2` | string | Yes | 3' base of the pair on the 3' side |
| `numStates` | integer | No | Equivalent states for degeneracy entropy (n=1 only, default 1) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `energy` | number | Bulge-loop free energy in kcal/mol (rounded to 2 dp) |

## Errors

| Code | Message |
|------|---------|
| 1001 | Expected a single character (length-1 string) |

## Examples

### Example 1: Single-nucleotide bulge with stacking

**User Prompt:**
> Free energy of a 1-nt bulge 'A' between G-C and G-C pairs.

**Expected Tool Call:**
```json
{
  "tool": "bulge_loop_energy",
  "arguments": {
    "bulgeSize": 1, "bulgedBase": "A",
    "pair5_base1": "G", "pair5_base2": "C",
    "pair3_base1": "G", "pair3_base2": "C"
  }
}
```

**Response:**
```json
{ "energy": 0.54 }
```
`init(1)=3.8 + stacking(GG/CC)=−3.26 = 0.54`.

### Example 2: Multi-nt bulge terminal penalty

**User Prompt:**
> 3-nt bulge between G-C and A-U pairs.

**Expected Tool Call:**
```json
{
  "tool": "bulge_loop_energy",
  "arguments": {
    "bulgeSize": 3, "bulgedBase": "A",
    "pair5_base1": "G", "pair5_base2": "C",
    "pair3_base1": "A", "pair3_base2": "U"
  }
}
```

**Response:**
```json
{ "energy": 3.65 }
```
`init(3)=3.2 + AU terminal penalty 0.45 = 3.65` (only the A-U closing pair is penalised).

## Performance

- **Time Complexity:** O(1).
- **Space Complexity:** O(1).

## See Also

- [internal_loop_energy](internal_loop_energy.md)
- [hairpin_loop_energy](hairpin_loop_energy.md)
- [stem_energy](stem_energy.md)
